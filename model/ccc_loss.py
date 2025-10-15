"""
Concordance Correlation Coefficient (CCC) Loss Module

This module implements the CCC loss function, an alternative to MSE for regression tasks.
CCC measures both precision and accuracy of predictions, making it suitable for
assessing agreement between predicted and true values.

The CCC ranges from -1 to 1:
    1: Perfect concordance
    0: No concordance
   -1: Perfect discordance

CCC combines:
    - Pearson correlation (precision/consistency)
    - Closeness to the line of perfect concordance (accuracy/bias)

Reference:
    Lin, L. I. (1989). "A concordance correlation coefficient to evaluate reproducibility."
    Biometrics, 45(1), 255-268.

Classes:
    CCCLoss: PyTorch loss module implementing CCC loss
"""

import torch


class CCCLoss(torch.nn.Module):
    """
    Concordance Correlation Coefficient (CCC) Loss Function.

    This loss is computed as 1 - CCC, so minimizing this loss maximizes concordance.

    The CCC is defined as:
        CCC = (2 * ρ * σ_x * σ_y) / (σ_x² + σ_y² + (μ_x - μ_y)²)

    where:
        ρ = Pearson correlation coefficient
        σ_x, σ_y = standard deviations of true and predicted values
        μ_x, μ_y = means of true and predicted values

    Components:
        Numerator: 2 * covariance
        Denominator: sum of variances + squared difference of means

    The CCC penalizes both:
        1. Poor correlation (lack of precision)
        2. Systematic bias (shift from identity line)

    Attributes:
        eps (float): Small constant for numerical stability (default: 1e-6)
    """

    def __init__(self, eps=1e-6):
        """
        Initialize the CCC loss module.

        Args:
            eps (float): Small epsilon value added for numerical stability
                        to prevent division by zero
        """
        super(CCCLoss, self).__init__()
        self.eps = eps

    def forward(self, y_true, y_hat):
        """
        Calculate the CCC loss between true and predicted values.

        Args:
            y_true (torch.Tensor): True values (ground truth labels)
            y_hat (torch.Tensor): Predicted values (model outputs)

        Returns:
            torch.Tensor: Scalar loss value (1 - CCC)
                         Lower values indicate better concordance
                         Range: [0, 2], optimal value = 0

        Example:
            >>> loss_fn = CCCLoss()
            >>> y_true = torch.tensor([1.0, 2.0, 3.0, 4.0])
            >>> y_pred = torch.tensor([1.1, 2.2, 2.9, 4.1])
            >>> loss = loss_fn(y_true, y_pred)
            >>> # loss will be close to 0 for good predictions
        """
        # Calculate means
        y_true_mean = torch.mean(y_true)
        y_hat_mean = torch.mean(y_hat)

        # Calculate variances
        y_true_var = torch.var(y_true)
        y_hat_var = torch.var(y_hat)

        # Calculate standard deviations
        y_true_std = torch.std(y_true)
        y_hat_std = torch.std(y_hat)

        # Center the data (remove means)
        vx = y_true - y_true_mean
        vy = y_hat - y_hat_mean

        # Calculate Pearson correlation coefficient (PCC)
        # PCC = covariance / (std_x * std_y)
        # Using numerically stable formula with epsilon
        pcc = torch.sum(vx * vy) / (
            torch.sqrt(torch.sum(vx ** 2) + self.eps) *
            torch.sqrt(torch.sum(vy ** 2) + self.eps)
        )

        # Calculate CCC
        # CCC = (2 * PCC * std_x * std_y) / (var_x + var_y + (mean_x - mean_y)²)
        ccc = (2 * pcc * y_true_std * y_hat_std) / (
            y_true_var + y_hat_var + (y_hat_mean - y_true_mean) ** 2
        )

        # Return loss as 1 - CCC
        # Perfect concordance (CCC = 1) gives loss = 0
        # No concordance (CCC = 0) gives loss = 1
        ccc_loss = 1 - ccc

        return ccc_loss
