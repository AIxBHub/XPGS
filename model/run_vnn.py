import argparse
import copy
import os
from data_wrapper import TrainingDataWrapper
from vnn_train import VNNTrainer

def main():
    #torch.set_printoptions(precision = 5)
    parser = argparse.ArgumentParser(description = 'Train VNN')
    parser.add_argument('-onto', help = 'Ontology file used to guide the neural network', type = str)
    parser.add_argument('-optimize', help = 'Hyper-parameter optimization (1=No optimization, 2=Optuna)', type = int, default = 1)
    parser.add_argument('-gene2id', help = 'Gene to ID mapping file', type = str)
    parser.add_argument('-min_dropout_layer', help = 'Start dropout from this Layer number', type = int, default = 2)
    parser.add_argument('-dropout_fraction', help = 'Dropout Fraction', type = float, default = 0.1)
    parser.add_argument('-lr', help = 'Learning rate', type = float, default = 0.0002)
    parser.add_argument('-wd', help = 'Weight decay', type = float, default =  0.00005)
    parser.add_argument('-epoch', help = 'Training epochs for training', type = int, default = 50)
    parser.add_argument('-alpha', help = 'Loss parameter alpha', type = float, default = 0.1)
    parser.add_argument('-batchsize', help = 'Batchsize', type = int, default = 64)
    parser.add_argument('-modeldir', help = 'Folder for trained models', type = str, default = 'MODEL/')
    parser.add_argument('-delta', help = 'Minimum change in loss to be considered an improvement', type = float, default = 0.01)
    parser.add_argument('-cuda', help = 'Specify GPU', type = int, default = 0)
    #parser.add_argument('-patience', help = 'Early stopping epoch limit', type = int, default = 30) # Used for optimization
    #parser.add_argument('-delta', help = 'Minimum change in loss to be considered an improvement', type = float, default = 0.001) # Used for optimization
    parser.add_argument('-train', help = 'Training dataset', type=str)
    parser.add_argument('-genotype_hiddens', help = 'Mapping for the number of neurons in each term in genotype parts', type = int, default = 50)
    parser.add_argument('-testsetratio', help = 'Assign test set ratio from training dataset if no test set given', type = float, default = 0.2)
    parser.add_argument('-test', help = 'Test dataset', type = str)
    parser.add_argument('-metric_output', help = 'Output table with loss and metrics every epoch', type=str, default="metrics_output.tsv")
    opt = parser.parse_args()
    data_wrapper = TrainingDataWrapper(opt)
    
    if not os.path.exists(opt.modeldir):
        os.makedirs(opt.modeldir)
    
    #if opt.optimize == 1:
    #    VNNTrainer(data_wrapper).train_model()

    #elif opt.optimize == 2:
    #    trial_params = OptunaNNTrainer(data_wrapper).exec_study()
    #    for key, value in trial_params.items():
    #        if hasattr(data_wrapper, key):
    #            setattr(data_wrapper, key, value)
    #    VNNTrainer(data_wrapper).train_model()

    #else:
    #    print("Wrong value for optimize.")
    #    exit(1)
    
if __name__ == "__main__":
	main()