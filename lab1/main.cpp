#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm> // Для std::round

const double PI = 3.14159265358979323846;

// Функция для генерации сигнала
std::vector<double> generateSignal(int numHarmonics, const std::vector<double>& amplitudes, double samplingRate, double initialPhase, int numPoints, double deltaT) {
    std::vector<double> signal(numPoints, 0.0);
    for (int i = 0; i < numPoints; ++i) {
        double t = i * deltaT;
        for (int h = 0; h < numHarmonics; ++h) {
            int harmonicFrequency = 50 * (h + 1); // Гармоники кратные 50 Гц
            signal[i] += amplitudes[h] * sin(2 * PI * harmonicFrequency * t + initialPhase);
        }
    }
    return signal;
}

// Функция для преобразования дискретного сигнала в цифровой (квантование)
std::vector<double> convertToDigital(const std::vector<double>& signal, int levels) {
    double maxAmplitude = *std::max_element(signal.begin(), signal.end());
    double minAmplitude = *std::min_element(signal.begin(), signal.end());
    double range = maxAmplitude - minAmplitude;

    std::vector<double> digitalSignal(signal.size());
    for (size_t i = 0; i < signal.size(); ++i) {
        // Квантование и масштабирование
        digitalSignal[i] = minAmplitude + (std::round((signal[i] - minAmplitude) / range * (levels - 1))) * (range / (levels - 1));
    }
    return digitalSignal;
}

// Функция для экспорта данных в файл
void exportToFile(const std::vector<double>& signal, const std::string& filename) {
    std::ofstream outFile(filename);
    for (double value : signal) {
        outFile << value << std::endl;
    }
    outFile.close();
}

int main() {
    // Начальные параметры для проверки
    int numHarmonics = 3; // Количество гармоник
    std::vector<double> amplitudes = {1.0, 0.5, 0.3}; // Амплитуды гармоник
    double samplingRate = 1000.0; // Частота дискретизации (Гц)
    double initialPhase = 0.0; // Начальный сдвиг (радианы)
    int numPoints = 100; // Количество точек (уменьшено для наглядности)
    double deltaT = 1.0 / samplingRate; // Шаг времени (секунды)
    int quantizationLevels = 16; // Уровни квантования для цифрового сигнала

    /*
    // Закомментированный блок для ввода параметров с клавиатуры
    int numHarmonics;
    std::cout << "Введите количество гармоник: ";
    std::cin >> numHarmonics;

    std::vector<double> amplitudes(numHarmonics);
    for (int i = 0; i < numHarmonics; ++i) {
        std::cout << "Введите амплитуду " << i + 1 << "-й гармоники: ";
        std::cin >> amplitudes[i];
    }

    double samplingRate;
    std::cout << "Введите частоту дискретизации: ";
    std::cin >> samplingRate;

    double initialPhase;
    std::cout << "Введите начальный сдвиг: ";
    std::cin >> initialPhase;

    int numPoints;
    std::cout << "Введите количество точек: ";
    std::cin >> numPoints;

    double deltaT;
    std::cout << "Введите шаг времени (дельта t): ";
    std::cin >> deltaT;

    int quantizationLevels;
    std::cout << "Введите количество уровней квантования: ";
    std::cin >> quantizationLevels;
    */

    // Генерация сигнала
    std::vector<double> analogSignal = generateSignal(numHarmonics, amplitudes, samplingRate, initialPhase, numPoints, deltaT);

    // Преобразование в цифровой сигнал
    std::vector<double> digitalSignal = convertToDigital(analogSignal, quantizationLevels);

    // Экспорт данных в файлы
    exportToFile(analogSignal, "analog_signal.txt");
    exportToFile(digitalSignal, "digital_signal.txt");

    std::cout << "Дискретный и цифровой сигналы успешно сгенерированы и экспортированы в файлы 'analog_signal.txt' и 'digital_signal.txt'" << std::endl;

    return 0;
}
