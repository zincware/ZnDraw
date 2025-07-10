import znsocket


def dict_to_nested_znsocket(data: dict, key: str, client: znsocket.Client, **kwargs) -> znsocket.Dict:
    dct = znsocket.Dict(client, key=key, **kwargs)
    for k, v in data.items():
        if isinstance(v, dict):
            dct[k] = dict_to_nested_znsocket(v, f"{key}.{k}", client, **kwargs)
        else:
            dct[k] = v
    return dct
