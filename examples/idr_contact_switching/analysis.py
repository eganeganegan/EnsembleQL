from pathlib import Path

import ensembleql as eql

HERE = Path(__file__).parent
trajectory = eql.load(HERE / "switching.xyz", topology=HERE / "switching.pdb")

switches = trajectory.query((HERE / "query.eql").read_text())
first_contact = trajectory.query("FIND CONTACT(resid 17, resid 42);")
second_contact = trajectory.query("FIND CONTACT(resid 17, resid 53);")

print("switches")
print(switches)
print("R17-D42 recurrence", first_contact.recurrence_statistics())
print("R17-E53 recurrence", second_contact.recurrence_statistics())
print("switch graph", switches.event_graph(within_ps=3000))
