from pathlib import Path

import ensembleql as eql

HERE = Path(__file__).parent
trajectory = eql.load(HERE / "switching.xyz", topology=HERE / "switching.pdb")
events = trajectory.query((HERE / "query.eql").read_text())
print(events)
