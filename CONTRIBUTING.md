# Contributing

Contributions that improve portability, test coverage, documentation, or bug fixes are welcome.

## Guidelines
- **Preserve scientific behavior:** Do not alter viral classification thresholds, database filtering rules, or taxonomy reconciliation logic without an explicit issue or pull request documenting the biological rationale.
- **No personal paths or credentials:** Never commit cluster-specific paths, email addresses, API keys, raw sequencing reads, databases, or compiled documents.
- **Small, focused pull requests:** Keep refactoring, dependency updates, and behavioral fixes in separate commits.

## Local Development Checks
Run the test and lint checks before opening a pull request:

```bash
# Fast stdlib unit, syntax, and hygiene checks (no external dependencies required)
python3 -m unittest discover -s tests -v

# Full suite including pandas/biopython tests and Ruff linting
pytest -v
ruff check tests workflow/phase1/scripts workflow/phase2/scripts

# Bash syntax check on all workflow shell scripts
for s in $(find workflow -name '*.sh' | sort); do bash -n "$s"; done
```
