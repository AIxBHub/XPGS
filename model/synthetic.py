import torch
import numpy as np

### GPT-4o generated ###

def generate_synthetic_data(num_samples=100, num_features=10, num_classes=2):
    # Generate random features to simulate gene expression data
    X = np.random.rand(num_samples, num_features).astype(np.float32)  # Simulate gene expression data with values between 0 and 1
    
    # Generate random binary labels for classification
    y = np.random.randint(num_classes, size=(num_samples,)).astype(np.int64)
    
    return torch.tensor(X), torch.tensor(y)

## Example usage
#X, y = generate_synthetic_data(num_samples=100, num_features=10)
#print("Features:", X[:5])  # Show the first 5 samples
#print("Labels:", y[:5])    # Show the corresponding labels