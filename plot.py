import plot_base

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.tri as tri

filename = "ref.msh"

[Nodes, Element] = plot_base.read(filename)

print(Nodes)
