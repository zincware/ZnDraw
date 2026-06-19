# ZnDraw JOSS Paper Design

**Date:** 2026-04-02
**Authors:** Fabian Zills, Rokas Elijošius
**Target:** Journal of Open Source Software (JOSS)
**Word limit:** 750–1750 words (target ~1100)

---

## Paper Structure

### 1. Summary (~180 words)

ZnDraw is a web-based, real-time collaborative visualization and editing tool for atomistic simulations. Built on a Python-first architecture, ZnDraw exposes trajectories as a `MutableSequence[ase.Atoms]` — any mutation, whether from Python or the built-in browser editor for moving atoms, is instantly synchronized to all connected clients via Socket.IO. This bidirectional sync enables researchers to inspect, manipulate, and analyze molecular trajectories both programmatically (scripts, notebooks) and interactively (3D browser GUI) simultaneously.

The frontend, built with React and Three.js, is shipped within a single cross-platform Python wheel — `pip install zndraw && zndraw file.xyz` requires no external services. ZnDraw supports 50+ file formats through ASE, lazy-streams large HDF5/H5MD trajectories, and renders 100,000+ atoms via GPU-instanced meshes. Interactive 2D Plotly figures complement the 3D scene for property visualization and analysis.

A Pydantic-based extension system allows users to register custom modifiers, selections, and analysis tools — with a set of custom frontend forms (e.g., Ketcher molecule editor, auto-generated property selectors) — that execute server-side through a FIFO worker queue. Extensions and data are scoped to rooms, enabling isolated collaborative sessions. For production, ZnDraw scales horizontally behind a load balancer with shared Redis state.

### 2. Statement of Need (~200 words)

The rise of Machine-Learned Interatomic Potentials (MLIPs) has made large-scale atomistic simulations accessible to broader research communities. As system sizes grow and collaborative workflows become common, researchers need tools that go beyond static desktop viewers — they need to visualize, edit, and analyze trajectories interactively while sharing results with collaborators in real time.

Existing molecular visualization tools force a fundamental trade-off. Desktop applications like OVITO, VMD, PyMOL, and ChimeraX offer rich functionality and scripting interfaces, but lack web accessibility and real-time collaboration. Web-based viewers such as nglview, py3Dmol, Chemiscope, and Mol* run in browsers but are read-only — users can view structures but cannot modify them, run custom analysis, or share a live session with collaborators.

No existing open-source tool combines real-time multi-user collaboration, a Python-native mutable API, and a server-side extension system in a web-based interface. ZnDraw addresses this gap. A researcher can load a trajectory from a Jupyter notebook, a collaborator can join the same room from a browser across the world, and a custom MLIP-based modifier can run server-side — all operating on the same synchronized state. ZnDraw is designed for the collaborative, Python-centric workflows that MLIP-driven research demands.

### 3. State of the Field (~200 words)

Molecular visualization has a long history of desktop tools. VMD uses a custom Tcl scripting language, PyMOL and ChimeraX offer Python scripting atop C/C++ cores, and OVITO provides Python bindings wrapping its C++ pipeline. While powerful, all are inherently single-user and single-machine.

Web-based alternatives have emerged but remain limited to viewing. nglview and py3Dmol wrap JavaScript libraries as Jupyter widgets for read-only display. Chemiscope provides interactive structure-property exploration but no editing or extensibility. Mol* powers the RCSB PDB viewer with high-performance biomolecular rendering but offers no Python API or collaboration. ASE's built-in X3D viewer provides minimal Jupyter display without trajectory playback or interaction.

Augmented and virtual reality tools have explored collaboration: chARpack enables co-located multi-user AR sessions on HoloLens hardware, and MolecularWebXR offers shared VR rooms via WebXR. However, both target immersive hardware rather than everyday computational workflows and lack Python APIs or integration with the atomistic simulation ecosystem.

ZnDraw is, to our knowledge, the only open-source tool that combines real-time multi-user collaboration via WebSockets, a Python-first mutable API built on ASE, a web-native browser interface, and a server-side extension system. This combination required a fundamentally different architecture — a networked, event-driven system with shared state — that no desktop viewer or read-only web widget could provide through incremental extension.

### 4. Software Design (~250 words)

ZnDraw's architecture follows a core design philosophy: every connected client — whether a web browser or Python client — must operate on the same synchronized state.

The backend is built on FastAPI with Socket.IO for bidirectional event-driven communication. Frame data — arbitrary per-atom and per-structure NumPy arrays — is serialized via MessagePack with NumPy support. Storage is layered: frame data is persisted through asebytes backends — in-memory by default, LMDB for single-host persistence, or MongoDB for distributed access. A SQL database (SQLite or PostgreSQL) manages users, rooms, and extension metadata. Redis handles locks, queues, and Socket.IO coordination. In standalone mode, all of these run in-process via fakeredis and in-memory SQLite, requiring no external services.

The Python client exposes trajectories as a `MutableSequence[ase.Atoms]`, allowing easy mutation of structures using familiar ASE-based workflows. Appending, slicing, or modifying frames from Python triggers Socket.IO events that update all connected browsers in real time, and vice versa for edits made in the browser's interactive atom editor. All state — frames, geometries, selections, figures — is scoped to rooms, providing isolation between collaborative sessions.

The extension system uses Pydantic models as the interface: users subclass `Extension`, declare parameters as typed fields, and implement a `run(vis)` method. The frontend renders custom forms (e.g., Ketcher for molecule input, auto-generated property selectors) based on the extension's schema. Extensions execute server-side through a FIFO Redis-backed worker queue, allowing long-running computations (e.g., MLIP relaxations) without blocking the UI. Multiple workers can consume from the queue for parallel throughput.

The React/Three.js frontend renders atoms via GPU-instanced meshes, enabling visualization of systems with 100,000+ atoms. It ships as static assets inside the Python wheel, making `pip install zndraw` a complete, cross-platform installation. For production, ZnDraw scales horizontally — multiple server instances behind a load balancer share state through Redis, with workers scaled independently.

### 5. Features and Implementation (~150 words + 2 code snippets + 1 figure)

**Figure:** Composite screenshot of the ZnDraw browser UI showing (1) 3D atom scene with the water/ethanol box, (2) Plotly figure panel, (3) extension sidebar with a form visible.

**Code snippet 1 — Python client with molify:**

```python
from zndraw import ZnDraw
from molify import smiles2conformers, pack

water = smiles2conformers("O", numConfs=2)
ethanol = smiles2conformers("CCO", numConfs=5)
box = pack([water, ethanol], [7, 5], density=1000)

vis = ZnDraw()
vis.append(box)
vis.selection = vis.select(smarts="[OH]")  # select all OH groups
```

**Code snippet 2 — Custom extension:**

```python
from zndraw.extensions import Extension, Category
from pydantic import Field

class Relaxation(Extension):
    category = Category.MODIFIER
    fmax: float = Field(0.05, description="Force convergence")

    def run(self, vis, **kwargs):
        from ase.optimize import LBFGS
        atoms = vis[vis.step]
        atoms.calc = kwargs["calculator"]
        LBFGS(atoms).run(fmax=self.fmax)
        vis.append(atoms)

vis.register_job(Relaxation, run_kwargs={"calculator": calc})
```

**Surrounding text:** Describe the CLI entry point (`zndraw file.xyz`), format support via ASE (50+ formats), lazy HDF5/H5MD streaming for large trajectories, room scoping for collaborative sessions, and interactive Plotly figure integration for property visualization.

### 6. Related Software (~100 words)

The functionality of ZnDraw builds upon and integrates with the following packages:

- **ASE**: For representing atomic structures and interfacing with simulation engines and file formats.
- **molify**: For SMILES-based molecule generation, conformer creation, and substructure selection.
- **Three.js / React**: For GPU-accelerated 3D rendering in the browser via instanced meshes.
- **Socket.IO**: For bidirectional real-time communication between server and clients.

ZnDraw is a core component in the following software packages:

- **IPSuite**: For interactive inspection of MLIP training workflows.
- **mlipx**: For visualizing MLIP benchmark results across real-world test scenarios.

### 7. Acknowledgements

F.Z. acknowledges support by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) in the framework of the priority program SPP 2363, "Utilization and Development of Machine Learning for Molecular Applications — Molecular Machine Learning" Project No. 497249646. Further funding through the DFG under Germany's Excellence Strategy – EXC 2075 – 390740016 and the Stuttgart Center for Simulation Science (SimTech) was provided.

**Note:** Check if Rokas has separate funding to acknowledge.

### 8. AI Usage Disclosure

One sentence required by JOSS review criteria (Jan 2026 policy). To be filled in at writing time — e.g., "Generative AI tools were [not] used in the development of this software."

### 9. References (~15–20 citations)

Key references to include:

- ASE (Larsen et al., 2017)
- molify (Zills, 2025 — JOSS)
- OVITO (Stukowski, 2010)
- VMD (Humphrey et al., 1996)
- PyMOL (Schrödinger)
- ChimeraX (Pettersen et al., 2021)
- nglview / NGL (Rose et al., 2018)
- py3Dmol / 3Dmol.js (Quer & Rego, 2015)
- Chemiscope (Fraux et al.)
- Mol* (Sehnal et al., 2021)
- chARpack (Rau et al., 2024 — DOI: 10.1021/acs.jcim.4c00462)
- MolecularWebXR
- Socket.IO
- FastAPI (Ramírez)
- Three.js
- React
- Redis
- PACKMOL (Martínez et al., 2009)
- RDKit (Landrum et al.)
- Elijošius et al. 2025 (NatComm — ZnDraw used for similarity kernels)
- IPSuite (Zills et al., 2024)

---

## Key Differentiators to Emphasize

1. **Only open-source tool with real-time multi-user collaboration** for atomistic visualization
2. **Python-first `MutableSequence[ase.Atoms]` API** — no other tool has this
3. **Pydantic extension system** with custom frontend forms, FIFO worker queue, multi-worker support
4. **Single `pip install`** ships React/Three.js frontend in the wheel — zero external services for standalone
5. **Horizontal scaling** — multiple instances behind load balancer, shared Redis state
6. **50+ formats via ASE**, lazy HDF5/H5MD streaming, 100k+ atom rendering via instancing
7. **Room-scoped isolation** for collaborative sessions

## JOSS Review Readiness Checklist

- [ ] License: EPL v2.0 (OSI-approved) — present in repo
- [ ] Tests: 1000+ tests, CI via GitHub Actions
- [ ] Documentation: README, API docs, developer guide
- [ ] Community guidelines: Need to verify CONTRIBUTING.md and CODE_OF_CONDUCT.md exist
- [ ] Installation: `pip install zndraw` works
- [ ] Development history: 6+ months of public commits
- [ ] Releases: Tagged versions via hatch-vcs

## Notes

- The NatComm paper (Elijošius et al., 2025) covers similarity kernels, not ZnDraw itself — no overlap with this JOSS submission
- JOSS does not require a formal "Research Impact Statement" section — impact evidence is evaluated during review, not as a paper section
- The molify JOSS paper can serve as a formatting template
