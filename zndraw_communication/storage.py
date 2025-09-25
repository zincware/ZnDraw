import numpy as np
import typing as t
import zarr.dtype
import json
from collections.abc import MutableSequence


def encode_data(data: dict) -> dict:
    serialized = {}
    for key, value in data.items():
        if isinstance(value, np.ndarray):
            serialized[key] = {
                'data': value.tobytes(),
                'shape': value.shape,
                'dtype': str(value.dtype)
            }
        elif isinstance(value, dict):
            serialized[key] = encode_data(value)
        else:
            serialized[key] = value
    return serialized

def decode_data(data: dict) -> dict:
    deserialized = {}
    for key, value in data.items():
        if isinstance(value, dict) and 'data' in value and 'shape' in value and 'dtype' in value:
            deserialized[key] = np.frombuffer(value['data'], dtype=value['dtype']).reshape(value['shape'])
        elif isinstance(value, dict):
            deserialized[key] = decode_data(value)
        else:
            deserialized[key] = value
    return deserialized

def create_zarr(root: zarr.Group, data: dict):
    for key, value in data.items():
        if isinstance(value, np.ndarray):
            initial_shape = (1,) + value.shape
            chunks = (1,) + value.shape
            dataset = root.create_array(
                name=key,
                shape=initial_shape,
                chunks=chunks,
                dtype=value.dtype
            )
            dataset[0] = value
        elif isinstance(value, dict):
            try:
                # if it is purely standard types,
                # we serialize it directly (assuming it being small)
                value = json.dumps(value)
                value = np.array(value)
                initial_shape = (1,) + value.shape
                chunks = (1,) + value.shape
                dataset = root.create_array(
                    name=key,
                    shape=initial_shape,
                    chunks=chunks,
                    dtype=zarr.dtype.VariableLengthUTF8()
                )
                dataset[0] = value
                dataset.attrs['format'] = "json"
            except TypeError:
                subgroup = root.create_group(key)
                create_zarr(subgroup, value)
        else:
            value = np.array(json.dumps(value))
            initial_shape = (1,) + value.shape
            chunks = (1,) + value.shape
            dataset = root.create_array(
                name=key,
                shape=initial_shape,
                chunks=chunks,
                dtype=zarr.dtype.VariableLengthUTF8()
            )
            dataset[0] = value
            dataset.attrs['format'] = "json"

def read_zarr(root: zarr.Group, index: int = 0) -> dict:
    data = {}
    for key in root:
        item = root[key]
        if isinstance(item, zarr.Array):
            if item.attrs.get('format') == "json":
                data[key] = json.loads(item[index].item())
            else:
                data[key] = item[index]
        elif isinstance(item, zarr.Group):
            data[key] = read_zarr(item, index=index)
    return data


def extend_zarr(root: zarr.Group, data: list[dict]):
    """
    Extends a Zarr hierarchy using the native zarr.Array.append() method.
    """
    def _append_recursive(group: zarr.Group, data_dict: dict):
        for key, value in data_dict.items():
            if key in group:
                item = group[key]
                if isinstance(item, zarr.Array):
                    if item.attrs.get('format') == 'json':
                        prepared_value = [np.array(json.dumps(value))]
                    else:
                        prepared_value = np.expand_dims(value, axis=0)
                    item.append(prepared_value)
                elif isinstance(item, zarr.Group):
                    _append_recursive(item, value)
            else:
                raise ValueError(f"Key '{key}' not found in Zarr group.")

    for entry in data:
        _append_recursive(root, entry)


class ZarrStorageSequence(MutableSequence):
    def __init__(self, group: zarr.Group):
        self.group = group

    def __getitem__(self, index: int|list[int]|slice) -> dict|list[dict]:
        if isinstance(index, slice):
            index = list(range(*index.indices(len(self))))
        is_single = False
        if isinstance(index, int):
            is_single = True
            index = [index]
        result = [read_zarr(self.group, i) for i in index]
        if is_single:
            return result[0]
        return result

    def __setitem__(self, index: int|list[int]|slice, value: dict|list[dict]):
        if isinstance(index, slice):
            index = list(range(*index.indices(len(self))))
        if isinstance(index, int):
            index = [index]

        raise NotImplementedError

    def __delitem__(self, index: int|list[int]|slice):
        if isinstance(index, slice):
            index = list(range(*index.indices(len(self))))
        if isinstance(index, int):
            index = [index]
        raise NotImplementedError

    def __len__(self):
        return self.group[list(self.group.keys())[0]].shape[0]

    def insert(self, index: int, value: dict):
        value = encode_data(value)
        raise NotImplementedError
    
    def append(self, value: dict) -> None:
        self.extend([value])

    def extend(self, values: list[dict]) -> None:
        extend_zarr(self.group, values)
