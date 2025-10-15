import argparse
import torch
import torch.nn as nn
import numpy as np
import pandas as pd
from data_wrapper import TrainingDataWrapper
import prepare_dataloader
import utils

class RLIPPCalculator:
    """
    Calculate RLIPP (Relative Local Improvement in Predictive Power) scores
    for each term/subsystem in the hierarchical neural network.

    RLIPP measures how much predictive power a term adds beyond its children.
    Positive RLIPP = term improves prediction relative to children
    Negative RLIPP = term decreases predictive power
    """

    def __init__(self, model, data_wrapper, device, predictor_epochs=100):
        self.model = model
        self.model.eval()
        self.data_wrapper = data_wrapper
        self.device = device
        self.predictor_epochs = predictor_epochs
        self.model.to(device)

    def get_embeddings_and_labels(self, data_loader):
        """
        Run forward pass on data and collect hidden embeddings for all terms
        and ground truth labels.
        """
        all_embeddings = {}
        all_labels = []

        with torch.no_grad():
            for i, (inputdata, labels) in enumerate(data_loader):
                print(f"Processing batch {i+1}/{len(data_loader)}")

                features = inputdata.to(torch.float32).to(self.device)
                labels = labels.to(torch.float32).unsqueeze(1).to(self.device)

                # Get all hidden embeddings from the model
                aux_out_map, hidden_embeddings_map = self.model(features)

                # Store embeddings for each term
                for term, embedding in hidden_embeddings_map.items():
                    if term not in all_embeddings:
                        all_embeddings[term] = []
                    all_embeddings[term].append(embedding.cpu())

                all_labels.append(labels.cpu())

        # Concatenate all batches
        for term in all_embeddings:
            all_embeddings[term] = torch.cat(all_embeddings[term], dim=0)

        all_labels = torch.cat(all_labels, dim=0)

        return all_embeddings, all_labels

    def predict_from_embedding(self, embedding):
        """
        Create prediction from an embedding using a simple linear layer.
        Returns predictions and the trained predictor.
        """
        # Simple linear predictor
        predictor = nn.Linear(embedding.shape[1], 1).to(self.device)
        optimizer = torch.optim.Adam(predictor.parameters(), lr=0.001)
        criterion = nn.MSELoss()

        return predictor

    def calculate_predictive_power(self, embedding, labels, epochs=100):
        """
        Calculate predictive power (correlation) of an embedding for the labels.
        Trains a simple linear model and returns the correlation.
        """
        embedding = embedding.to(self.device)
        labels = labels.to(self.device)

        # Create simple linear predictor
        predictor = nn.Linear(embedding.shape[1], 1).to(self.device)
        optimizer = torch.optim.Adam(predictor.parameters(), lr=0.001)
        criterion = nn.MSELoss()

        # Quick training
        predictor.train()
        for epoch in range(epochs):
            optimizer.zero_grad()
            predictions = predictor(embedding)
            loss = criterion(predictions, labels)
            loss.backward()
            optimizer.step()

        # Evaluate
        predictor.eval()
        with torch.no_grad():
            predictions = predictor(embedding)
            corr = utils.pearson_corr(predictions, labels)

        return corr.item()

    def calculate_rlipp_scores(self, data_loader):
        """
        Calculate RLIPP score for each term in the hierarchy.

        For each term:
        1. Calculate predictive power using the term's embedding
        2. Calculate predictive power using only children's embeddings (concatenated)
        3. RLIPP = (power_with_term - power_with_children) / max(abs(power_with_term), abs(power_with_children))
        """
        print("Extracting embeddings from model...")
        all_embeddings, all_labels = self.get_embeddings_and_labels(data_loader)

        rlipp_scores = {}

        # Iterate through each layer from bottom to top
        print("\nCalculating RLIPP scores for each term...")
        for layer_idx, layer in enumerate(self.model.term_layer_list):
            print(f"\nProcessing layer {layer_idx + 1}/{len(self.model.term_layer_list)}")

            for term_idx, term in enumerate(layer):
                print(f"  Term {term_idx + 1}/{len(layer)}: {term}")

                # Get embedding for current term
                term_embedding = all_embeddings[term]

                # Get children embeddings
                children = self.model.term_neighbor_map[term]

                if len(children) == 0:
                    # Leaf node - no children to compare against
                    # Calculate predictive power using term embedding only
                    term_power = self.calculate_predictive_power(term_embedding, all_labels, epochs=self.predictor_epochs)
                    rlipp_scores[term] = {
                        'rlipp': term_power,  # For leaf nodes, RLIPP is just the predictive power
                        'term_power': term_power,
                        'children_power': 0.0,
                        'num_children': 0,
                        'layer': layer_idx
                    }
                else:
                    # Concatenate children embeddings
                    children_embeddings = [all_embeddings[child] for child in children]
                    combined_children = torch.cat(children_embeddings, dim=1)

                    # Calculate predictive power with term embedding
                    term_power = self.calculate_predictive_power(term_embedding, all_labels, epochs=self.predictor_epochs)

                    # Calculate predictive power with only children embeddings
                    children_power = self.calculate_predictive_power(combined_children, all_labels, epochs=self.predictor_epochs)

                    # Calculate RLIPP (relative improvement)
                    # Normalize by the maximum absolute value to get relative improvement
                    max_power = max(abs(term_power), abs(children_power))
                    if max_power > 0:
                        rlipp = (term_power - children_power) / max_power
                    else:
                        rlipp = 0.0

                    rlipp_scores[term] = {
                        'rlipp': rlipp,
                        'term_power': term_power,
                        'children_power': children_power,
                        'num_children': len(children),
                        'layer': layer_idx
                    }

                    print(f"    Term power: {term_power:.4f}, Children power: {children_power:.4f}, RLIPP: {rlipp:.4f}")

        return rlipp_scores

    def save_rlipp_scores(self, rlipp_scores, output_file):
        """
        Save RLIPP scores to a TSV file, sorted by RLIPP score.
        """
        # Convert to dataframe
        rows = []
        for term, scores in rlipp_scores.items():
            rows.append({
                'term': term,
                'rlipp': scores['rlipp'],
                'term_power': scores['term_power'],
                'children_power': scores['children_power'],
                'num_children': scores['num_children'],
                'layer': scores['layer']
            })

        df = pd.DataFrame(rows)
        df = df.sort_values('rlipp', ascending=False)
        df.to_csv(output_file, sep='\t', index=False)

        print(f"\nRLIPP scores saved to: {output_file}")
        print(f"\nTop 10 terms by RLIPP score:")
        print(df.head(10).to_string(index=False))

        print(f"\nBottom 10 terms by RLIPP score:")
        print(df.tail(10).to_string(index=False))

        return df


def main():
    parser = argparse.ArgumentParser(description='Calculate RLIPP scores from trained VNN model')
    parser.add_argument('-model', help='Path to saved model (.pt file)', type=str, required=True)
    parser.add_argument('-onto', help='Ontology file used to guide the neural network', type=str, required=True)
    parser.add_argument('-gene2id', help='Gene to ID mapping file', type=str, required=True)
    parser.add_argument('-test', help='Test dataset (.pt file)', type=str, required=True)
    parser.add_argument('-batchsize', help='Batch size', type=int, default=64)
    parser.add_argument('-output', help='Output file for RLIPP scores', type=str, default='rlipp_scores.tsv')
    parser.add_argument('-cuda', help='Specify GPU', type=int, default=0)
    parser.add_argument('-genotype_hiddens', help='Number of neurons in each term (must match training)', type=int, default=50)
    parser.add_argument('-min_dropout_layer', help='Min dropout layer (must match training)', type=int, default=2)
    parser.add_argument('-dropout_fraction', help='Dropout fraction (must match training)', type=float, default=0.1)
    parser.add_argument('-predictor_epochs', help='Training epochs for predictive power calculation', type=int, default=100)

    args = parser.parse_args()

    # Set device
    device = torch.device(f"cuda:{args.cuda}" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # Load the trained model
    print(f"\nLoading model from: {args.model}")
    model = torch.load(args.model, map_location=device)
    model.eval()

    print("\nModel loaded successfully!")
    print(f"Model has {len(model.term_layer_list)} layers")
    print(f"Root term: {model.root}")

    # Create a minimal data wrapper for loading test data
    # We need this to maintain compatibility with the model structure
    class MinimalArgs:
        pass

    minimal_args = MinimalArgs()
    minimal_args.gene2id = args.gene2id
    minimal_args.onto = args.onto
    minimal_args.genotype_hiddens = args.genotype_hiddens
    minimal_args.min_dropout_layer = args.min_dropout_layer
    minimal_args.dropout_fraction = args.dropout_fraction
    minimal_args.lr = 0.001  # Not used
    minimal_args.wd = 0.0001  # Not used
    minimal_args.alpha = 0.1  # Not used
    minimal_args.epoch = 1  # Not used
    minimal_args.batchsize = args.batchsize
    minimal_args.modeldir = './'  # Not used
    minimal_args.delta = 0.01  # Not used
    minimal_args.train = args.test  # Use test as train (we just need the structure)
    minimal_args.test = args.test
    minimal_args.testsetratio = 0.0
    minimal_args.optimize = 1
    minimal_args.metric_output = 'temp.tsv'

    data_wrapper = TrainingDataWrapper(minimal_args, logger=None)

    # Load test data
    print(f"\nLoading test data from: {args.test}")
    _, test_loader = prepare_dataloader.get_from_pt(args.test, args.test, args.batchsize)

    # Calculate RLIPP scores
    calculator = RLIPPCalculator(model, data_wrapper, device, predictor_epochs=args.predictor_epochs)

    rlipp_scores = calculator.calculate_rlipp_scores(test_loader)

    # Save results
    df = calculator.save_rlipp_scores(rlipp_scores, args.output)

    # Print summary statistics
    print("\n" + "="*60)
    print("RLIPP Score Summary Statistics:")
    print("="*60)
    print(df['rlipp'].describe())
    print(f"\nTerms with positive RLIPP (improve prediction): {(df['rlipp'] > 0).sum()}")
    print(f"Terms with negative RLIPP (decrease prediction): {(df['rlipp'] < 0).sum()}")
    print(f"Terms with near-zero RLIPP (|RLIPP| < 0.01): {(df['rlipp'].abs() < 0.01).sum()}")


if __name__ == "__main__":
    main()
