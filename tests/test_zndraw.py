from zndraw import io


def test_io_get_graph(water):
    G = io.get_graph(water)

    assert G.number_of_nodes() == 3
    assert G.nodes[1] == {
        "symbol": "H",
        "number": 1,
        "x": 0.0,
        "y": 0.763239,
        "z": -0.477047,
    }
    assert G.nodes[2] == {
        "symbol": "H",
        "number": 1,
        "x": 0.0,
        "y": -0.763239,
        "z": -0.477047,
    }
    assert G.nodes[0] == {"symbol": "O", "number": 8, "x": 0.0, "y": 0.0, "z": 0.119262}

    assert len(G.edges) == 2
    assert G.has_edge(0, 1)
    assert G.has_edge(0, 2)
    assert not G.has_edge(1, 2)
