### Objective
# The objective of this script is to generate plots for results from goes_sensible_heat_flux.py

### Imports
import matplotlib.pyplot as plt

### Main functions
def main(iter_list, L_list, err_list):
    plotter(iter_list, L_list, err_list)

if __name__ == "__main__":
    main()
    
def plotter(iter_list, L_list, err_list):
    # Plot can be updated with variables for each field for easier editing/replacement        
    plt.plot(iter_list[2:], L_list[2:])
    plt.plot(iter_list[2:], err_list[2:])
    plt.title('Convergence statistics')
    plt.xlabel('Iteration #')
    plt.ylabel('L (m), Error (%)')
    plt.show()