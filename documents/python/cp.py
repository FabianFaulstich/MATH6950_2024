import numpy as np
from tensorly.decomposition import parafac

# Create a tensor
X = np.array([[[1, 2], [3, 4]],
              [[5, 6], [7, 8]]])

print("Original Tensor:")
print(X)

# Perform CP decomposition with rank 1
factors_rank1 = parafac(X, rank=1)
A1, B1, C1 = factors_rank1.factors

# Perform CP decomposition with rank 2
factors_rank2 = parafac(X, rank=2)
A2, B2, C2 = factors_rank2.factors

# Reconstruct the tensor from the factors
X_reconstructed_rank1 = np.einsum('ia, jb, kc -> ijk', A1, B1, C1)
X_reconstructed_rank2 = np.einsum('ia, jb, kc -> ijk', A2, B2, C2)

print("\nFactors for CP decomposition with rank 1:")
print("Factor A:")
print(A1)
print("Factor B:")
print(B1)
print("Factor C:")
print(C1)
print("\nReconstructed Tensor (Rank 1):")
print(X_reconstructed_rank1)

print("\nFactors for CP decomposition with rank 2:")
print("Factor A:")
print(A2)
print("Factor B:")
iprint(B2)
print("Factor C:")
print(C2)
print("\nReconstructed Tensor (Rank 2):")
print(X_reconstructed_rank2)

