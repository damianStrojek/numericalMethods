import csv
from matplotlib import pyplot as plt

nodes_numbers = [4, 5, 6, 8, 12, 21, 31]
#nodes_numbers = [2, 3, 4, 7, 8, 11, 21]
methods = ["Lagrange", "splines"]
routes = {
    "ostrowa": "Gdańsk Ostrowa",
    "genoa_rapallo": "Genoa-Rapallo",
    "diff_heights": "Different heights",
    "tczew_starogard": "Tczew Starogard"
}

for filename, route in routes.items():
    for method in methods:
        for nodes_number in nodes_numbers:
            # interpolation nodes
            x_nodes = []
            y_nodes = []
            with open('C:\\Users\\Administrator\\source\\repos\\profilWysokosciowy'
                      '\\results\\{}\\{}_{}_nodes.txt'.format(method, filename, nodes_number), 'r') as file:
                reader = csv.reader(file, delimiter=' ')
                for row in reader:
                    x_nodes.append(float(row[0]))
                    y_nodes.append(float(row[1]))

            # results of interpolation
            x = []
            y = []
            with open('C:\\Users\\Administrator\\source\\repos\\profilWysokosciowy'
                      '\\results\\{}\\{}_{}.txt'.format(method, filename, nodes_number), 'r') as file:
                reader = csv.reader(file, delimiter=' ')
                for row in reader:
                    x.append(float(row[0]))
                    y.append(float(row[1]))

            # reference values
            x_ref = []
            y_ref = []
            with open('C:\\Users\\Administrator\\source\\repos\\profilWysokosciowy'
                      '\\data\\{}.txt'.format(filename), 'r') as file:
                reader = csv.reader(file, delimiter=' ')
                for row in reader:
                    if float(row[0]) <= x[len(x) - 1]:
                        x_ref.append(float(row[0]))
                        y_ref.append(float(row[1]))

            fig = plt.figure()
            plt.title("{} – {} method\n{} nodes".format(route, method, nodes_number))
            plt.xlabel("distance [m]")
            plt.ylabel("elevation [m]")
            plt.plot(x_ref, y_ref, color='blue', label='actual values')
            plt.plot(x, y, color='red', label='interpolated values')
            plt.plot(x_nodes, y_nodes, 'o', color='red')
            plt.legend(loc='best', borderaxespad=0.)
            plt.savefig('C:\\Users\\Administrator\\source\\repos\\profilWysokosciowy'
                        '\\charts\\{}\\{}_{}_chart.png'.format(method, filename, nodes_number))
            plt.close(fig)