from pathlib import Path

import ensembleql as eql

HERE = Path(__file__).parent
trajectory = eql.load(HERE / "adsorption.xyz", topology=HERE / "adsorption.pdb")
query = (HERE / "query.eql").read_text()
events = trajectory.query(query)

print(events)
print(events.recurrence_statistics())
