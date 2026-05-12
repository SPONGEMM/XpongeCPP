# Reference Xponge Force-Field Snapshot

This directory is a read-only data snapshot copied from the reference Xponge
repository path `Xponge/forcefield`.

Amber runtime data now lives under `src/XpongeCPP/data/amber` as the installable
package resource used by the active C++ workflow.

This `third_party/xponge_reference_forcefield` tree is not the primary runtime
path. It exists as a source-tree reference snapshot for non-Amber migration,
baseline comparison, and parser bring-up without depending on a separate
checkout of the original Xponge repository.
