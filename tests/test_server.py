import time

import ase
import ase.collections
import pytest
from ase.build import molecule

from zndraw import ZnDraw


@pytest.mark.usefixtures("setup")
class TestZnDraw:
    def test_title(self, server):
        self.driver.get(server)
        assert self.driver.title == "ZnDraw"

    def test_vis_connection(self, server):
        self.driver.get(server)
        vis = ZnDraw(url=server)
        assert vis.socket.connected

    def test_vis_len(self, server):
        self.driver.get(server)
        time.sleep(1)
        # we need to wait for all the data to be loaded.
        # This includes jsonschemas and atoms.
        vis = ZnDraw(url=server)
        assert len(vis) == 1
        assert vis[0] == ase.Atoms()

        vis[0] = molecule("H2O")
        assert len(vis) == 1
        assert vis[0] == molecule("H2O")

    # def test_vis_selection(self, server):
    # TODO: fix
    #     self.driver.get(server)
    #     time.sleep(1)
    #     vis = ZnDraw(url=server)
    #     vis[0] = molecule("H2O")
    #     vis.selection = [1, 2]
    #     assert vis.selection == [1, 2]

    def test_insert_before(self, server):
        self.driver.get(server)
        time.sleep(1)
        s22 = list(ase.collections.s22)
        vis = ZnDraw(url=server)
        vis.extend(s22)
        del vis[0]
        assert len(vis) == 22
        assert list(vis) == s22

        vis.step = 10
        assert len(vis) == 22
        vis.insert(0, s22[5])
        assert len(vis) == 23
        # assert vis.step == 11 # TODO: fix vis.step assertion
        assert vis[0] == s22[5]
        assert vis[1] == s22[0]
        assert vis[2] == s22[1]

    def test_insert_after(self, server):
        self.driver.get(server)
        time.sleep(1)
        s22 = list(ase.collections.s22)
        vis = ZnDraw(url=server)
        vis.extend(s22)
        del vis[0]
        assert len(vis) == 22
        assert list(vis) == s22

        vis.step = 10
        assert len(vis) == 22
        vis.insert(11, s22[5])
        assert len(vis) == 23
        # assert vis.step == 10
        assert vis[10] == s22[10]
        assert vis[11] == s22[5]
        assert vis[12] == s22[11]

    def test_append_to_non_empty(self, server):
        self.driver.get(server)
        time.sleep(1)
        vis = ZnDraw(url=server)
        vis.append(molecule("H2O"))
        vis.append(molecule("H2O"))
        assert len(vis) == 3
        assert vis[0] == ase.Atoms()
        assert vis[1] == molecule("H2O")
        assert vis[2] == molecule("H2O")

    def test_extend_to_non_empty(self, server):
        self.driver.get(server)
        time.sleep(1)
        vis = ZnDraw(url=server)
        vis.extend([molecule("H2O"), molecule("H2O")])
        assert len(vis) == 3
        assert vis[0] == ase.Atoms()
        assert vis[1] == molecule("H2O")
        assert vis[2] == molecule("H2O")

    def test_delete_backwards(self, server):
        self.driver.get(server)
        time.sleep(1)
        s22 = list(ase.collections.s22)
        vis = ZnDraw(url=server)
        # vis[:] = list(ase.collections.s22) # not supported
        # vis.extend(ase.collections.s22) # not working, because first element is already there
        for idx, atoms in enumerate(s22):
            vis[idx] = atoms

        assert len(vis) == 22
        for idx in range(22):
            assert vis[idx] == s22[idx]

        for idx in range(22):
            del vis[len(vis) - 1]
            assert len(vis) == 22 - idx - 1
            for jdx in range(len(vis)):
                assert vis[jdx] == s22[jdx]

    def test_delete_forwards(self, server):
        self.driver.get(server)
        time.sleep(1)
        s22 = list(ase.collections.s22)
        vis = ZnDraw(url=server)
        # vis[:] = list(ase.collections.s22) # not supported
        # vis.extend(ase.collections.s22) # not working, because first element is already there
        for idx, atoms in enumerate(s22):
            vis[idx] = atoms

        assert len(vis) == 22
        for idx in range(22):
            assert vis[idx] == s22[idx]

        for idx in range(22):
            del vis[0]
            assert len(vis) == 22 - idx - 1
            for jdx in range(len(vis)):
                assert vis[jdx] == s22[jdx + idx + 1]

    def test_delete_middle(self, server):
        self.driver.get(server)
        time.sleep(1)
        s22 = list(ase.collections.s22)
        vis = ZnDraw(url=server)
        # vis[:] = list(ase.collections.s22) # not supported
        # vis.extend(ase.collections.s22) # not working, because first element is already there
        for idx, atoms in enumerate(s22):
            vis[idx] = atoms

        assert len(vis) == 22
        for idx in range(22):
            assert vis[idx] == s22[idx]

        for idx in range(0, 22):
            if len(vis) <= 10:
                with pytest.raises(IndexError):
                    del vis[10]
            else:
                del vis[10]
                assert len(vis) == 22 - idx - 1
                for jdx in range(len(vis)):
                    assert vis[jdx] == s22[jdx] if jdx < 10 else s22[jdx + 1]

    def test_delete_slice_backwards(self, server):
        self.driver.get(server)
        time.sleep(1)
        s22 = list(ase.collections.s22)
        vis = ZnDraw(url=server)
        # vis[:] = list(ase.collections.s22) # not supported
        # vis.extend(ase.collections.s22) # not working, because first element is already there
        for idx, atoms in enumerate(s22):
            vis[idx] = atoms

        assert len(vis) == 22
        for idx in range(22):
            assert vis[idx] == s22[idx]

        for idx in range(22):
            del vis[len(vis) - 1 :]
            assert len(vis) == 22 - idx - 1
            for jdx in range(len(vis)):
                assert vis[jdx] == s22[jdx]

    def test_delete_slice_forwards(self, server):
        self.driver.get(server)
        time.sleep(1)
        s22 = list(ase.collections.s22)
        vis = ZnDraw(url=server)
        # vis[:] = list(ase.collections.s22) # not supported
        # vis.extend(ase.collections.s22) # not working, because first element is already there
        for idx, atoms in enumerate(s22):
            vis[idx] = atoms

        assert len(vis) == 22
        for idx in range(22):
            assert vis[idx] == s22[idx]

        for idx in range(22):
            del vis[:1]
            assert len(vis) == 22 - idx - 1
            for jdx in range(len(vis)):
                assert vis[jdx] == s22[jdx + idx + 1]

    def test_delete_slice_middle(self, server):
        self.driver.get(server)
        time.sleep(1)
        s22 = list(ase.collections.s22)
        vis = ZnDraw(url=server)
        # vis[:] = list(ase.collections.s22) # not supported
        # vis.extend(ase.collections.s22) # not working, because first element is already there
        for idx, atoms in enumerate(s22):
            vis[idx] = atoms

        assert len(vis) == 22
        for idx in range(22):
            assert vis[idx] == s22[idx]

        for idx in range(0, 22):
            del vis[10:11]
            if len(vis) <= 10:
                assert len(vis) == 10
                for jdx in range(len(vis)):
                    assert vis[jdx] == s22[jdx] if jdx < 10 else s22[jdx + 1]
            else:
                assert len(vis) == 22 - idx - 1
                for jdx in range(len(vis)):
                    assert vis[jdx] == s22[jdx] if jdx < 10 else s22[jdx + 1]

    def test_delete_all(self, server):
        self.driver.get(server)
        time.sleep(1)
        s22 = list(ase.collections.s22)
        vis = ZnDraw(url=server)
        # vis[:] = list(ase.collections.s22) # not supported
        # vis.extend(ase.collections.s22) # not working, because first element is already there
        for idx, atoms in enumerate(s22):
            vis[idx] = atoms

        assert len(vis) == 22
        for idx in range(22):
            assert vis[idx] == s22[idx]

        del vis[:]
        assert len(vis) == 0
        assert vis[:] == []
