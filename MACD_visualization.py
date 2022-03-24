# Damian Strojek s184407
import csv
from math import floor
from matplotlib import pyplot
from dateutil import parser


class Macd:
    __macd_series = []
    __signal_series = []
    __buy_sell_signals = []
    __buy_sell_signals_delayed = []
    __samples = []

    def __init__(self, n, samples):
        self.__samples = samples
        self.calculate_macd(n, samples)
        self.calculate_signal()
        self.calculate_buy_sell_signals()
        self.calculate_buy_sell_signals_delayed()

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

    def calculate_buy_sell_signals(self):
        self.__buy_sell_signals = ["nothing"]
        tmp_macd = self.get_macd(35)
        tmp_signals = self.get_signal()
        for i in range(1, len(tmp_macd)):
            if tmp_macd[i - 1] > tmp_signals[i - 1] and tmp_macd[i] <= tmp_signals[i]:
                self.__buy_sell_signals.append("sell")
            elif tmp_macd[i - 1] < tmp_signals[i - 1] and tmp_macd[i] >= tmp_signals[i]:
                self.__buy_sell_signals.append("buy")
            else:
                self.__buy_sell_signals.append("nothing")

    def calculate_buy_sell_signals_delayed(self):
        self.__buy_sell_signals_delayed = ["nothing"]
        tmp_macdD = self.get_macd(35)
        tmp_signalsD = self.get_signal()
        for i in range(2, len(tmp_macdD)):
            if tmp_macdD[i - 2] > tmp_signalsD[i] and tmp_macdD[i] <= tmp_signalsD[i]:
                self.__buy_sell_signals_delayed.append("sell")
            elif tmp_macdD[i - 2] < tmp_signalsD[i] and tmp_macdD[i] >= tmp_signalsD[i]:
                self.__buy_sell_signals_delayed.append("buy")
            else:
                self.__buy_sell_signals_delayed.append("nothing")

    def get_macd(self, i=0):
        return self.__macd_series[i::]

    def get_signal(self, i=0):
        return self.__signal_series[i::]

    def get_buy_sell_signals(self, i=0):
        return self.__buy_sell_signals[i::]

    def get_buy_sell_signals_delayed(self, i=0):
        return self.__buy_sell_signals_delayed[i::]

    def debug(self):
        for i in self.__samples:
            print(i)


def diagram_macd(time, macd, signal):
    pyplot.plot(time, macd, label="MACD", color="green")
    pyplot.plot(time, signal, label="Signal", color="red")
    pyplot.legend()
    pyplot.grid(True)
    pyplot.ylabel('Wartość składowych')
    pyplot.xlabel('Data (yyyy-mm)')
    pyplot.title('Wskaźnik MACD')
    pyplot.show()

# not working as I would like it to work
def diagram_macd_buy_sell(time, macd, signal, buy_sell):
    pyplot.plot(time, macd, label="MACD", color="green")
    pyplot.plot(time, signal, label="Signal", color="red")
    pyplot.legend()
    pyplot.grid(True)
    pyplot.ylabel('Wartość składowych')
    pyplot.xlabel('Data (yyyy-mm)')
    pyplot.title('Wskaźnik MACD BUY SELL')
    for i in range(1, len(buy_sell)):
        x = [time[i], -5]
        y = [time[i], 1]
        if buy_sell[i] == "buy":
            pyplot.plot(x, y, 'b', linewidth=1, color="blue")
        elif buy_sell[i] == "sell":
            pyplot.plot(x, y, 'b', linewidth=1, color="blue")
    pyplot.show()


def diagram_samples(time, data):
    pyplot.plot(time, data, label="macd", color='green')
    pyplot.ylabel('Kurs USD')
    pyplot.xlabel('Data (yyyy-mm)')
    pyplot.title('Kurs akcji Coca-Cola')
    pyplot.grid(True)
    pyplot.show()


def read_data(filename):
    csvfile = open(filename, newline='')
    cocacola_reader = list(csv.reader(csvfile, delimiter=','))

    # order : date, opening price, highest price, lowest price,
    # closing price, adjusted closing price, volume
    time = [parser.parse(i[0]) for i in cocacola_reader]
    samples = [float(j[4]) for j in cocacola_reader]
    return time, samples


def simulation(samples, signals, money):
    currency = 0
    for i in range(len(signals)):
        if signals[i] == "buy":
            if money != 0:
                currency = float(money / samples[i])
                money = 0
        elif signals[i] == "sell":
            if currency != 0:
                money = currency * samples[i]
                currency = 0
    if money == 0:
        # -1 is the last element
        return float(currency * samples[-1])
    else:
        return money


if __name__ == '__main__':
    n = 1000
    starting_money = 1000
    time, samples = read_data('cocacola.csv')
    macd = Macd(n, samples)

    print("Kapitał początkowy : " + str(starting_money))
    money_simulation = simulation(samples[35::], macd.get_buy_sell_signals(), starting_money)
    profit = floor(money_simulation / starting_money * 100) - 100
    print("Kapitał końcowy z algorytmu podstawowego : " + str(money_simulation))
    print("Profit : " + str(profit) + "%")

    print("Kapitał początkowy : " + str(starting_money))
    money_simulation = simulation(samples[35::], macd.get_buy_sell_signals_delayed(), starting_money)
    profit = floor(money_simulation / starting_money * 100) - 100
    print("Kapitał końcowy z opóźnionego algorytmu MACD : " + str(money_simulation))
    print("Profit : " + str(profit) + "%")

    diagram_macd(time[35::], macd.get_macd(35), macd.get_signal())
    diagram_samples(time[35::], samples[35::])
    #diagram_macd_buy_sell(time[35::], macd.get_macd(35), macd.get_signal(), macd.get_buy_sell_signals())

