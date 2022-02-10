#!/usr/bin/env python
# coding=utf-8
#!/usr/bin/env python
# coding=utf-8
"""
Model architecture
"""
import torch
import torch.nn as nn
import torch.nn.functional as F
from torchvision import transforms

from models import SingleImageBreastModel


class ResNet22(nn.Module):

    def __init__(self, L=256, D=128, k=4, batch_size=16, dropout=True, pretrained=True):
        super(ResNet22, self).__init__()
        self.k = k
        self.single_model = SingleImageBreastModel(1)
        self.single_model.load_state_from_shared_weights(
                state_dict=torch.load("/breast_cancer_classifier/models/ImageOnly__ModeImage_weights.p")["model"],
                view='L-CC',
                )
        self.global_avg_pool = nn.AdaptiveAvgPool2d(1)
        self.feature_extractor = self.single_model.view_resnet
       
        self.classifier = nn.Linear(1024, 3)
        self.batch_size = batch_size

    def forward(self, x):
    
        x = x.view(-1, 1, 512, 512)
        x = self.feature_extractor(x)
        x = self.global_avg_pool(x)
        x = x.view(x.shape[0], x.shape[1])
        x = x.view(-1, 4, 256)
        x = x.view(-1, 4*256)
        out = self.classifier(x)
        return out


class BMGatedAttention(nn.Module):

    def __init__(self):
        super(BMGatedAttention, self).__init__()
        self.D = 128
        self.K = 1

        feature_extractor = resnet34(True)
        feature_extractor.conv1 = nn.Conv2d(1, 64, kernel_size=(7,7),
                                            stride=(2,2), padding=(3,3),
                                            bias=False)
        self.feature_extractor = nn.Sequential(*(list(feature_extractor.children())[:-1]))
        self.attention_V = nn.Sequential(nn.Linear(512, self.D), nn.Tanh())
        self.attention_U = nn.Sequential(nn.Linear(512, self.D), nn.Sigmoid())
        self.attention_weights = nn.Linear(self.D, self.K)
        self.classifier = nn.Sequential(nn.Linear(512*self.K, 1), nn.Sigmoid())

    def forward(self, x):
        x = x.view(-1, 1, 512, 512)

        H = self.feature_extractor(x)
        H = H.reshape(-1, 512)
        H_ = H.view(-1, 4, 512)
        H_ = H_.permute(0, 2, 1)

        A_V = self.attention_V(H)
        A_U = self.attention_U(H)
        A = self.attention_weights(A_V * A_U)
        A = A.view(-1, 4, 1)
        A = F.softmax(A, dim=1)

        M = torch.bmm(H_, A)
        M = M.squeeze(-1)

        Y_prob = self.classifier(M)
        Y_hat = torch.ge(Y_prob, 0.5).float()
        return Y_prob, Y_hat, A


if __name__ == "__main__":
   
