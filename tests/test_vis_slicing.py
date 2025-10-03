from zndraw.zndraw import ZnDraw
import pytest
import numpy as np

def test_getitem_slice(server, s22):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    # Basic slicing
    slice1 = vis[5:15]
    assert len(slice1) == 10
    for i in range(10):
        assert slice1[i] == s22[5 + i]
    slice2 = vis[:10]
    assert len(slice2) == 10
    for i in range(10):
        assert slice2[i] == s22[i]
    slice3 = vis[10:]
    assert len(slice3) == len(s22) - 10
    for i in range(len(s22) - 10):
        assert slice3[i] == s22[10 + i]
    slice4 = vis[-10:]
    assert len(slice4) == 10
    for i in range(10):
        assert slice4[i] == s22[len(s22) - 10 + i]
    slice5 = vis[:-10]
    assert len(slice5) == len(s22) - 10
    for i in range(len(s22) - 10):
        assert slice5[i] == s22[i]
    
    # Slicing with step
    slice6 = vis[0:20:2]
    assert len(slice6) == 10
    for i in range(10):
        assert slice6[i] == s22[0 + i * 2]
    slice7 = vis[1:20:3]
    assert len(slice7) == 7
    for i in range(7):
        assert slice7[i] == s22[1 + i * 3]
    slice8 = vis[::5]
    assert len(slice8) == (len(s22) + 4) // 5
    for i in range(len(slice8)):
        assert slice8[i] == s22[i * 5]
    slice9 = vis[::-1]
    assert len(slice9) == len(s22)
    for i in range(len(s22)):
        assert slice9[i] == s22[len(s22) - 1 - i]
    slice10 = vis[-1:-11:-1]
    assert len(slice10) == 10
    for i in range(10):
        assert slice10[i] == s22[len(s22) - 1 - i]

def test_getitem_slice_invalid(server, s22):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    # Invalid slice types
    with pytest.raises(TypeError):
        _ = vis["invalid"]
    with pytest.raises(TypeError):
        _ = vis[5.5]
    with pytest.raises(TypeError):
        _ = vis[5:15:1.5]

    # Out of range slices should not raise errors
    slice1 = vis[1000:1010]
    assert len(slice1) == 0
    slice2 = vis[-1000:-990]
    assert len(slice2) == 0

    # Slices with step 0 should raise ValueError
    with pytest.raises(ValueError):
        _ = vis[::0]
    with pytest.raises(ValueError):
        _ = vis[5:15:0]

def test_getitem_index(server, s22):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    # Valid indices
    for i in range(len(s22)):
        assert vis[i] == s22[i]
    for i in range(-len(s22), 0):
        assert vis[i] == s22[len(s22) + i]
    for i in range(len(s22)):
        assert vis[i] == s22[np.array(i)]

    # Invalid indices
    with pytest.raises(IndexError):
        _ = vis[len(s22)]
    with pytest.raises(IndexError):
        _ = vis[-len(s22) - 1]

def test_getitem_list(server, s22):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    slice = vis[[1, 2, 3]]
    assert len(slice) == 3
    for i in range(3):
        assert slice[i] == s22[i + 1]
    
    slice = vis[[0, -1, len(s22) - 2]]
    assert len(slice) == 3
    assert slice[0] == s22[0]
    assert slice[1] == s22[len(s22) - 1]
    assert slice[2] == s22[len(s22) - 2]

    slice = vis[np.array([0, 2, -1])]
    assert len(slice) == 3
    assert slice[0] == s22[0]
    assert slice[1] == s22[2]
    assert slice[2] == s22[len(s22) - 1]

    with pytest.raises(IndexError):
        _ = vis[[0, len(s22)]]
    with pytest.raises(IndexError):
        _ = vis[[-len(s22) - 1, 0]]
    with pytest.raises(TypeError):
        _ = vis[[1, 2.5, 3]]
    with pytest.raises(TypeError):
        _ = vis[[1, "2", 3]]
    with pytest.raises(TypeError):
        _ = vis[[1, [2], 3]]