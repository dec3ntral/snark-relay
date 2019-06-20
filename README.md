# snark-relay

PoC implementation of a Bitcoin relay for Ethereum using zkSNARKS for verification.

Currently, zkSNARKS are used to verify:
- Merkle tree inclusion proofs (~990k gates for a max tree depth of 11)
