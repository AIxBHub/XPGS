import numpy as np
import time
import torch
import torch.nn as nn
import torch.optim as optim
import torch.utils.data as du
from torch.autograd import Variable

import utils
import prepare_dataloader
from dcell_nn import DCellNN
from bb_nn import blackbox
from ccc_loss import *
import os

class VNNTrainer():

    def __init__(self, data_wrapper):
        self.data_wrapper = data_wrapper


    def train_model(self):

        self.model = DCellNN(self.data_wrapper)## model build with GO ontology processed in data_wrapper
#        self.model.cuda(self.data_wrapper.cuda)

        min_loss = None
        max_corr = None
        term_mask_map = utils.create_term_mask(self.model.term_direct_gene_map, self.model.gene_dim, self.data_wrapper.cuda)
        for name, param in self.model.named_parameters():
            term_name = name.split('_')[0]
            if '_direct_gene_layer.weight' in name:
                param.data = torch.mul(param.data, term_mask_map[term_name]) #* 0.1
            else:
                param.data = param.data #* 0.1

        train_loader, val_loader = prepare_dataloader.get_data(self.data_wrapper.train, self.data_wrapper.test, self.data_wrapper.testsetratio, self.data_wrapper.gene_id_mapping, self.data_wrapper.batchsize)

        optimizer = torch.optim.AdamW(self.model.parameters(), lr=self.data_wrapper.lr, betas=(0.9, 0.99), eps=1e-05, weight_decay=self.data_wrapper.wd)
        optimizer.zero_grad()

        #print("epoch\ttrain_corr\ttrain_loss\ttrue_auc\tpred_auc\tval_corr\tval_loss\tgrad_norm\telapsed_time")
        logout = open(self.data_wrapper.outfile, "w")
        logout.write("epoch\ttrain_corr\ttrain_loss\ttrue_auc\tpred_auc\tval_corr\tval_loss\telapsed_time\n")
        
        epoch_start_time = time.time()
        for epoch in range(self.data_wrapper.epochs):
            # Train
            self.model.train()
            train_predict = torch.zeros(0, 0).cuda(self.data_wrapper.cuda)
            _gradnorms = torch.empty(len(train_loader)).cuda(self.data_wrapper.cuda) # tensor for accumulating grad norms from each batch in this epoch

            for i, (inputdata, labels) in enumerate(train_loader):
                # Convert torch tensor to Variable
                features = utils.build_input_vector(inputdata, self.data_wrapper.gene_id_mapping)
                features = features.to(torch.float32).cuda(self.data_wrapper.cuda)
                labels = labels.to(torch.float32).unsqueeze(1).cuda(self.data_wrapper.cuda)

                # Forward + Backward + Optimize
                optimizer.zero_grad()  # zero the gradient buffer
                cur_out = self.model(features)
                aux_out_map, _ = cur_out

                if train_predict.size()[0] == 0:
                    train_predict = aux_out_map['final'].data
                    train_label_gpu = labels
                else:
                    train_predict = torch.cat([train_predict, aux_out_map['final'].data], dim=0)
                    train_label_gpu = torch.cat([train_label_gpu, labels], dim=0)

                total_loss = 0
                for name, output in aux_out_map.items():
                    #loss = CCCLoss()
                    loss = nn.MSELoss()
                    if name == 'final':
                        total_loss += loss(output, labels)
                    else:
                        total_loss += self.data_wrapper.alpha * loss(output, labels)

                total_loss.backward()

                for name, param in self.model.named_parameters():
                    if '_direct_gene_layer.weight' not in name:
                        continue
                    term_name = name.split('_')[0]
                    param.grad.data = torch.mul(param.grad.data, term_mask_map[term_name])

                #_gradnorms[i] = utils.get_grad_norm(self.model.parameters(), 2.0).unsqueeze(0) # Save gradnorm for batch
                optimizer.step()

            gradnorms = sum(_gradnorms).unsqueeze(0).cpu().numpy()[0] # Save total gradnorm for epoch
            train_corr = utils.pearson_corr(train_predict, train_label_gpu)

            self.model.eval()

            val_predict = torch.zeros(0, 0).cuda(self.data_wrapper.cuda)

            val_loss = 0
            for i, (inputdata, labels) in enumerate(val_loader):
                # Convert torch tensor to Variable
                vfeatures = utils.build_input_vector(inputdata, self.data_wrapper.gene_id_mapping)
                vfeatures = vfeatures.to(torch.float32).cuda(self.data_wrapper.cuda)
                labels = labels.to(torch.float32).unsqueeze(1).cuda(self.data_wrapper.cuda)

                aux_out_map, hidden_embeddings_map = self.model(vfeatures)

                if val_predict.size()[0] == 0:
                    val_predict = aux_out_map['final'].data
                    val_label_gpu = labels
                else:
                    val_predict = torch.cat([val_predict, aux_out_map['final'].data], dim=0)
                    val_label_gpu = torch.cat([val_label_gpu, labels], dim=0)

                for name, output in aux_out_map.items():
                    #loss = CCCLoss()
                    loss = nn.MSELoss()
                    if name == 'final':
                        val_loss += loss(output, labels)

            val_corr = utils.pearson_corr(val_predict, val_label_gpu)

            epoch_end_time = time.time()
            true_auc = torch.mean(train_label_gpu)
            pred_auc = torch.mean(train_predict)
            #print("{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}".format(epoch, train_corr, total_loss, true_auc, pred_auc, val_corr, val_loss, gradnorms, epoch_end_time - epoch_start_time))
            print(f"epoch{epoch}, val_corr={val_corr:3f}")
            logout.write("{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\n".format(epoch, train_corr, total_loss, true_auc, pred_auc, val_corr, val_loss, epoch_end_time - epoch_start_time))
            epoch_start_time = epoch_end_time


            """if min_loss == None:
                min_loss = val_loss
                torch.save(self.model, self.data_wrapper.modeldir + '/model_final.pt')
                print("Model saved at epoch {}".format(epoch))
            elif min_loss - val_loss > self.data_wrapper.delta:
                min_loss = val_loss
                torch.save(self.model, self.data_wrapper.modeldir + '/model_final.pt')
                print("Model saved at epoch {}".format(epoch))"""
            pt_path = os.path.join(self.data_wrapper.modeldir, f'{self.data_wrapper.strfime}_model.pt')
            
            if max_corr == None and val_corr > 0.9:
                max_corr = val_corr
                torch.save(self.model,pt_path)
                print("Model saved at epoch {}".format(epoch))
            elif max_corr and (val_corr-max_corr) > self.data_wrapper.delta:
                max_corr = val_corr
                torch.save(self.model,pt_path)
                print("Model saved at epoch {}".format(epoch))
                
        # Calculate rlipp
        print("get feature in each term")
        #get_feature    
        return min_loss

