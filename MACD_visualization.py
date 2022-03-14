import csv
import os
from matplotlib import pyplot
from dateutil import parser

MONEY = 1000

class Macd:
    __macd_series = []
    __signal_series = []
    __buy_sell_signals = []
    __samples = []

    def __init__(self, n, samples):
        self.__samples = samples
        self.calculate_macd(n, samples)
        self.calculate_signal()
        # self.calculate_buy_sell_signals()

    def calculate_ema(self, period, samples, day):
        alpha = 2 / (period + 1)
        temp = samples[day - period: day + 1:]
        temp.reverse()
        numerator = float(0.0)
        denominator = float(0.0)
        for i in range(period + 1):
            numerator += (1 - alpha) ** i * temp[i]
            denominator += (1 - alpha) ** i
        return numerator / denominator

    def calculate_macd(self, n, samples):
        self.__macd_series = []
        for i in range(n):
            if i >= 26:
                ema12 = self.calculate_ema(12, samples, i)
                ema26 = self.calculate_ema(26, samples, i)
                self.__macd_series.append(ema12 - ema26)
            else:
                self.__macd_series.append(float(0))

    def calculate_signal(self):
        self.__signal = []
        for i in range(35, len(self.__macd_series)):
            self.__signal_series.append(self.calculate_ema(9, self.__macd_series, i))

    def get_macd(self, i = 0):
        return self.__macd_series[i::]

    def get_signal(self, i = 0):
        return self.__signal_series[i::]

    def debug(self):
        for i in self.__samples:
            print(i)


class Visualization:
    def diagram_macd(self, time, macd, signal):
        pyplot.plot(time, macd, label="MACD", color="blue")
        pyplot.plot(time, signal, label="Signal", color="red")
        pyplot.legend()
        pyplot.grid(True)
        pyplot.ylabel('Wartość składowych')
        pyplot.xlabel('Data (yyyy-mm)')
        pyplot.title('Wskaźnik MACD')
        pyplot.show()


def read_data(filename):
    csvfile = open(filename, newline='')
    cocacola_reader = list(csv.reader(csvfile, delimiter=','))
    cocacola_reader.reverse()

    # order : date, opening price, highest price, lowest price,
    # closing price, adjusted closing price, volume
    time = [parser.parse(i[0]) for i in cocacola_reader]
    samples = [float(j[4]) for j in cocacola_reader]
    return time, samples


if __name__ == '__main__':
    n = 1000
    time, samples = read_data('MACD_cocacola.csv')
    macd = Macd(n, samples)
    Visualization.diagram_macd(Visualization(), time[35::], macd.get_macd(35), macd.get_signal())
