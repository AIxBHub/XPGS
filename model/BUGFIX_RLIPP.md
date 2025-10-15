# Bug Fix: calculate_rlipp.py AttributeError

## Issue

**Error Message:**
```
AttributeError: 'method' object has no attribute 'epochs'
```

**Location:** Line 264 in `calculate_rlipp.py`

## Root Cause

The code was trying to set an attribute on a method object:
```python
calculator.calculate_predictive_power.epochs = args.predictor_epochs  # WRONG!
```

This is incorrect because:
- `calculate_predictive_power` is a method, not an object with settable attributes
- The `epochs` parameter should be passed when calling the method, not set as an attribute

## Solution

### Changes Made

1. **Modified `RLIPPCalculator.__init__()` to accept `predictor_epochs` parameter:**
   ```python
   def __init__(self, model, data_wrapper, device, predictor_epochs=100):
       self.model = model
       self.model.eval()
       self.data_wrapper = data_wrapper
       self.device = device
       self.predictor_epochs = predictor_epochs  # Store as instance attribute
       self.model.to(device)
   ```

2. **Updated `calculate_rlipp_scores()` to use stored `predictor_epochs`:**
   ```python
   # For leaf nodes
   term_power = self.calculate_predictive_power(
       term_embedding,
       all_labels,
       epochs=self.predictor_epochs  # Pass as parameter
   )

   # For parent nodes
   term_power = self.calculate_predictive_power(
       term_embedding,
       all_labels,
       epochs=self.predictor_epochs
   )
   children_power = self.calculate_predictive_power(
       combined_children,
       all_labels,
       epochs=self.predictor_epochs
   )
   ```

3. **Fixed instantiation in `main()` function:**
   ```python
   # OLD (WRONG):
   calculator = RLIPPCalculator(model, data_wrapper, device)
   calculator.calculate_predictive_power.epochs = args.predictor_epochs

   # NEW (CORRECT):
   calculator = RLIPPCalculator(
       model,
       data_wrapper,
       device,
       predictor_epochs=args.predictor_epochs
   )
   ```

4. **Fixed `TrainingDataWrapper` initialization:**
   ```python
   # Added logger=None parameter to avoid issues with new logging system
   data_wrapper = TrainingDataWrapper(minimal_args, logger=None)
   ```

## Files Modified

- [calculate_rlipp.py](calculate_rlipp.py)
  - Line 20: Added `predictor_epochs` parameter to `__init__()`
  - Line 25: Store `predictor_epochs` as instance attribute
  - Lines 135, 149, 152: Pass `epochs=self.predictor_epochs` to method calls
  - Line 264: Fixed calculator instantiation
  - Line 257: Added `logger=None` to `TrainingDataWrapper` call

## Testing

To verify the fix works:

```bash
# Test with default predictor epochs (100)
python model/calculate_rlipp.py \
    -model MODEL/trained_model.pt \
    -onto data/ontology.txt \
    -gene2id data/gene2ind.txt \
    -test data/test.pt \
    -output rlipp_scores.tsv

# Test with custom predictor epochs
python model/calculate_rlipp.py \
    -model MODEL/trained_model.pt \
    -onto data/ontology.txt \
    -gene2id data/gene2ind.txt \
    -test data/test.pt \
    -output rlipp_scores.tsv \
    -predictor_epochs 50  # Custom value
```

## Verification

The fix ensures:
- ✅ No `AttributeError` when running the script
- ✅ The `predictor_epochs` argument is properly passed to all `calculate_predictive_power()` calls
- ✅ Users can customize the number of epochs for predictive power calculation
- ✅ Default value of 100 epochs is used if not specified

## Related Changes

This fix is compatible with the logging enhancements made to the training scripts, which is why the `logger=None` parameter was added to the `TrainingDataWrapper` initialization.
