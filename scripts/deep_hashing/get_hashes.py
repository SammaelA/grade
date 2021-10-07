import os
import sys
import argparse
import warnings
import numpy as np
import scipy.io as sio
import model.dch as model
import data_provider.image as dataset

from pprint import pprint

warnings.filterwarnings("ignore")

#os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

sys.argv = ["get_hashes"]
parser = argparse.ArgumentParser(description='Triplet Hashing')
parser.add_argument('--lr', '--learning-rate', default=0.005, type=float)
parser.add_argument('--output-dim', default=64, type=int)   # 256, 128
parser.add_argument('--alpha', default=0.5, type=float)
parser.add_argument('--bias', default=0.0, type=float)
parser.add_argument('--gamma', default=20, type=float)
parser.add_argument('--iter-num', default=2000, type=int)
parser.add_argument('--q-lambda', default=0, type=float)
parser.add_argument('--dataset', default='tmp', type=str)
parser.add_argument('--gpus', default='0', type=str)
parser.add_argument('--log-dir', default='tflog', type=str)
parser.add_argument('-b', '--batch-size', default=128, type=int)
parser.add_argument('-vb', '--val-batch-size', default=8, type=int)
parser.add_argument('--decay-step', default=10000, type=int)
parser.add_argument('--decay-factor', default=0.1, type=int)

tanh_parser = parser.add_mutually_exclusive_group(required=False)
tanh_parser.add_argument('--with-tanh', dest='with_tanh', action='store_true')
tanh_parser.add_argument('--without-tanh', dest='with_tanh', action='store_false')
parser.set_defaults(with_tanh=True)

parser.add_argument('--img-model', default='alexnet', type=str)
parser.add_argument('--model-weights', type=str,
                    default='./models/tm_1_dch.npy')
parser.add_argument('--finetune-all', default=True, type=bool)
parser.add_argument('--save-dir', default="./models/", type=str)
parser.add_argument('--data-dir', default="./", type=str)
parser.add_argument('-e', '--evaluate', dest='evaluate', action='store_true')

args = parser.parse_args([])

os.environ['CUDA_VISIBLE_DEVICES'] = args.gpus

label_dims = {'cifar10': 10, 'cub': 200, 'nuswide_81': 81, 'coco': 80, 'tmp': 33}
Rs = {'cifar10': 54000, 'nuswide_81': 5000, 'coco': 5000, 'tmp' : 1000}#13200
args.R = Rs[args.dataset]
args.label_dim = label_dims[args.dataset]

args.filename = os.path.join(args.data_dir, args.dataset, "prepare.txt")

pprint(vars(args))

data_root = os.path.join(args.data_dir, args.dataset)

images = dataset.Dataset('img', data_root, args.filename, train=False)
arr = model.get_hashes(images, args)

arr_size_x = arr.shape[0]
arr_size_y = arr.shape[1]
arr_bytes = arr.tobytes()