# Retention guard: outputs later analysis needs must not be removed.
# The protected set and the reasons are in configs/retained_outputs.tsv;
# the logic is scripts/retention_guard.py.
#
# Included LAST from Snakefile_phase2, so every rule is defined by the time the
# check runs. It runs at every parse -- dry runs included -- and stops the
# workflow before anything executes if a rule declares a retained output
# temp() or deletes one in its shell that it does not itself produce.
#
# onsuccess/onerror run the ledger audit, which reports any retained file that
# existed after an earlier run and has since disappeared (by whatever means:
# a rule, a script, a manual rm). The audit cannot fail a finished run, so it
# reports loudly instead; `python3 scripts/retention_guard.py audit` exits 1.
import subprocess as _subprocess

sys.path.insert(0, os.path.join(workflow.basedir, "scripts"))
from retention_guard import check_workflow as _retention_check

_retention_check(workflow)


def _retention_audit():
    try:
        _subprocess.run([sys.executable, "scripts/retention_guard.py", "audit"], check=False)
    except Exception as e:  # never mask the run's own result
        print(f"[retention] audit could not run: {e}", file=sys.stderr)


onsuccess:
    _retention_audit()

onerror:
    _retention_audit()
