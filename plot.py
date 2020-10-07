### Objective
# The objective of this script is to generate plots for results from goes_sensible_heat_flux.py

### Imports
import matplotlib.pyplot as plt
import numpy as np

### Main functions
def main(i, dataset_names, *args, **kwargs):
    plotter(i, dataset_names, *args, **kwargs)

if __name__ == "__main__":
    main()

### Plotter accepts:
# i (int): number of iterations on outer loop of goes_sensible_heat_flux.py
# dataset_names (list): list of names for the datasets plotted in args, kwargs
# *args: list items, typically from the val_dict, to plot values for parameters of interest
# **kwargs: list items, typically from the error_dict, to plot error values for parameters of interest    
def plotter(i, dataset_names, *args, **kwargs):
    # Initialize plot 
    fig, ax1 = plt.subplots()
    
    # Plot metadata
    plt.title('Model run statistics')
    plt.xlabel('Iteration #')
    
    # Iterand for legend label selection. Used due to the secondary y-axis
    n = 0
    ax1.tick_params(axis='y')  
    # Loop to plot all value data in *args
    for arg in args:      
        ax1.plot(np.linspace(0, i, i), arg, label=dataset_names[n])   
        n += 1
        
    ax2 = ax1.twinx()
    ax2.tick_params(axis='y')
    ax2.set_ylabel('Error')
    # Loop to plot all error value data in **kwargs
    for val in kwargs.values():
        ax2.plot(np.linspace(0, i, i), val, label=dataset_names[n])
        n += 1        
    
    # Re-group each y-axis label, combine into one legens
    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(lines, [line.get_label() for line in lines])
        
    plt.show()