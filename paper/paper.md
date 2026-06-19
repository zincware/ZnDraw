---
title: 'ZnDraw: Real-Time Collaborative Visualization for Atomistic Simulations'
tags:
  - Python
  - visualization
  - MLIPs
  - ASE
  - collaboration
  - Three.js
  - FastAPI
authors:
  - name: Fabian Zills
    orcid: 0000-0002-6936-4692
    affiliation: "1"
  - name: Rokas Elijošius
    orcid: 0000-0000-0000-0000
    affiliation: "2"
affiliations:
 - name: Institute for Computational Physics, University of Stuttgart, 70569 Stuttgart, Germany
   index: 1
 - name: TODO
   index: 2
date: 2 April 2026
bibliography: bibliography.bib
---

# Summary

ZnDraw is a web-based, real-time collaborative visualization and editing tool for atomistic simulations.
Built on a Python-first architecture, ZnDraw exposes trajectories as a `MutableSequence[ase.Atoms]` --- any mutation, whether from Python or the built-in browser editor for moving atoms, is instantly synchronized to all connected clients via Socket.IO.
This bidirectional sync enables researchers to inspect, manipulate, and analyze molecular trajectories both programmatically (scripts, notebooks) and interactively (3D browser GUI) simultaneously.

The frontend, built with React and Three.js, is shipped within a single cross-platform Python wheel --- `pip install zndraw && zndraw file.xyz` requires no external services.
ZnDraw supports over 50 file formats through ASE [@larsenAtomicSimulationEnvironment2017], lazy-streams large HDF5/H5MD trajectories, and renders 100,000+ atoms via GPU-instanced meshes.
Interactive 2D Plotly figures complement the 3D scene for property visualization and analysis.

A Pydantic-based extension system allows users to register custom modifiers, selections, and analysis tools --- with a set of custom frontend forms (e.g., Ketcher molecule editor, auto-generated property selectors) --- that execute server-side through a FIFO worker queue.
Extensions and data are scoped to rooms, enabling isolated collaborative sessions.
For production, ZnDraw scales horizontally behind a load balancer with shared Redis state.

# Statement of need

The rise of Machine-Learned Interatomic Potentials (MLIPs) has made large-scale atomistic simulations accessible to broader research communities.
As system sizes grow and collaborative workflows become common, researchers need tools that go beyond static desktop viewers --- they need to visualize, edit, and analyze trajectories interactively while sharing results with collaborators in real time.

Existing molecular visualization tools force a fundamental trade-off.
Desktop applications like OVITO [@stukowskiVisualizationAnalysisAtomistic2010], VMD [@humphreyVMDVisualMolecular1996], PyMOL [@schrodingerPyMOL], and ChimeraX [@pettersenUCSFChimeraXStructure2021] offer rich functionality and scripting interfaces, but lack web accessibility and real-time collaboration.
Web-based viewers such as nglview [@nguyenNGLviewInteractiveMolecular2018], py3Dmol [@querPy3DmolEnablingGlance2015], Chemiscope [@frauxChemiscope2020], and Mol\* [@sehnalMolStarStateArt2021] run in browsers but are read-only --- users can view structures but cannot modify them, run custom analysis, or share a live session with collaborators.

No existing open-source tool combines real-time multi-user collaboration, a Python-native mutable API, and a server-side extension system in a web-based interface.
ZnDraw addresses this gap.
A researcher can load a trajectory from a Jupyter notebook, a collaborator can join the same room from a browser across the world, and a custom MLIP-based modifier can run server-side --- all operating on the same synchronized state.
ZnDraw is designed for the collaborative, Python-centric workflows that MLIP-driven research demands.

# State of the field

Molecular visualization has a long history of desktop tools.
VMD [@humphreyVMDVisualMolecular1996] uses a custom Tcl scripting language, PyMOL [@schrodingerPyMOL] and ChimeraX [@pettersenUCSFChimeraXStructure2021] offer Python scripting atop C/C++ cores, and OVITO [@stukowskiVisualizationAnalysisAtomistic2010] provides Python bindings wrapping its C++ pipeline.
While powerful, all are inherently single-user and single-machine.

Web-based alternatives have emerged but remain limited to viewing.
nglview [@nguyenNGLviewInteractiveMolecular2018] and py3Dmol [@querPy3DmolEnablingGlance2015] wrap JavaScript libraries as Jupyter widgets for read-only display.
Chemiscope [@frauxChemiscope2020] provides interactive structure-property exploration but no editing or extensibility.
Mol\* [@sehnalMolStarStateArt2021] powers the RCSB PDB viewer with high-performance biomolecular rendering but offers no Python API or collaboration.
ASE's [@larsenAtomicSimulationEnvironment2017] built-in X3D viewer provides minimal Jupyter display without trajectory playback or interaction.

Augmented and virtual reality tools have explored collaboration: chARpack [@rauChARpackChemistryAugmented2024] enables co-located multi-user AR sessions on HoloLens hardware, and MolecularWebXR [@cassianoMolecularWebXR2025] offers shared VR rooms via WebXR.
However, both target immersive hardware rather than everyday computational workflows and lack Python APIs or integration with the atomistic simulation ecosystem.

ZnDraw is, to our knowledge, the only open-source tool that combines real-time multi-user collaboration via WebSockets, a Python-first mutable API built on ASE, a web-native browser interface, and a server-side extension system.
This combination required a fundamentally different architecture --- a networked, event-driven system with shared state --- that no desktop viewer or read-only web widget could provide through incremental extension.

# Software design

ZnDraw's architecture follows a core design philosophy: every connected client --- whether a web browser or Python client --- must operate on the same synchronized state.

The backend is built on FastAPI with Socket.IO for bidirectional event-driven communication.
Frame data --- arbitrary per-atom and per-structure NumPy arrays --- is serialized via MessagePack with NumPy support.
Storage is layered: frame data is persisted through asebytes backends --- in-memory by default, LMDB for single-host persistence, or MongoDB for distributed access.
A SQL database (SQLite or PostgreSQL) manages users, rooms, and extension metadata.
Redis handles locks, queues, and Socket.IO coordination.
In standalone mode, all of these run in-process via fakeredis and in-memory SQLite, requiring no external services.

The Python client exposes trajectories as a `MutableSequence[ase.Atoms]`, allowing easy mutation of structures using familiar ASE-based workflows.
Appending, slicing, or modifying frames from Python triggers Socket.IO events that update all connected browsers in real time, and vice versa for edits made in the browser's interactive atom editor.
All state --- frames, geometries, selections, figures --- is scoped to rooms, providing isolation between collaborative sessions.

The extension system uses Pydantic models as the interface: users subclass `Extension`, declare parameters as typed fields, and implement a `run(vis)` method.
The frontend renders custom forms (e.g., Ketcher for molecule input, auto-generated property selectors) based on the extension's schema.
Extensions execute server-side through a FIFO Redis-backed worker queue, allowing long-running computations (e.g., MLIP relaxations) without blocking the UI.
Multiple workers can consume from the queue for parallel throughput.

The React/Three.js frontend renders atoms via GPU-instanced meshes, enabling visualization of systems with 100,000+ atoms.
It ships as static assets inside the Python wheel, making `pip install zndraw` a complete, cross-platform installation.
For production, ZnDraw scales horizontally --- multiple server instances behind a load balancer share state through Redis, with workers scaled independently.

# Features and implementation

![The ZnDraw browser interface showing a 3D visualization of a water--ethanol mixture (left), an interactive Plotly figure for property analysis (right), and the extension sidebar (bottom).\label{fig:zndraw-ui}](zndraw_ui.png)

A typical ZnDraw workflow begins with creating molecular structures using `molify` [@zillsMolifyMolecularStructure2025] and visualizing them interactively:

```python
from zndraw import ZnDraw
from molify import smiles2conformers, pack

water = smiles2conformers("O", numConfs=2)
ethanol = smiles2conformers("CCO", numConfs=5)
box = pack([water, ethanol], [7, 5], density=1000)

vis = ZnDraw()
vis.append(box)
vis.selection = vis.select(smarts="[OH]")
```

Custom extensions integrate seamlessly into the UI and execute server-side:

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

ZnDraw can also be launched from the command line with `zndraw file.xyz`, supporting any of the over 50 file formats recognized by ASE.
For large trajectories stored in HDF5 or H5MD format, frames are lazy-streamed without loading the entire file into memory.
Interactive Plotly figures linked to per-frame properties provide complementary 2D analysis alongside the 3D scene, as shown in \autoref{fig:zndraw-ui}.

# Related software

The functionality of ZnDraw builds upon and integrates with the following packages:

- [ASE](https://wiki.fysik.dtu.dk/ase/) [@larsenAtomicSimulationEnvironment2017]: For representing atomic structures and interfacing with simulation engines and file formats.
- [molify](https://github.com/zincware/molify) [@zillsMolifyMolecularStructure2025]: For SMILES-based molecule generation, conformer creation, and substructure selection.
- [Three.js](https://threejs.org/): For GPU-accelerated 3D rendering in the browser via instanced meshes.
- [Socket.IO](https://socket.io/): For bidirectional real-time communication between server and clients.

ZnDraw is a core component in the following software packages:

- [IPSuite](https://github.com/zincware/ipsuite) [@zillsCollaborationMachineLearnedPotentials2024]: For interactive inspection of MLIP training workflows.
- [mlipx](https://github.com/basf/mlipx): For visualizing MLIP benchmark results across real-world test scenarios.

# Acknowledgements

F.Z. acknowledges support by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) in the framework of the priority program SPP 2363, "Utilization and Development of Machine Learning for Molecular Applications --- Molecular Machine Learning" Project No. 497249646. Further funding through the DFG under Germany's Excellence Strategy -- EXC 2075 -- 390740016 and the Stuttgart Center for Simulation Science (SimTech) was provided.

# References
