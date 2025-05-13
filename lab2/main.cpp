#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include <numeric>

const double PI = 3.14159265358979323846;
const double ERROR_THRESHOLD = 0.05; // 5% погрешность

// Функция для генерации сигнала
std::vector<double> generateSignal(int numHarmonics, const std::vector<double>& amplitudes,
                                  double samplingRate, double initialPhase, int numPoints, double deltaT) {
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
        digitalSignal[i] = minAmplitude + (std::round((signal[i] - minAmplitude) / range * (levels - 1))) * (range / (levels - 1));
    }
    return digitalSignal;
}

// Прямое преобразование Фурье (DFT)
std::vector<std::complex<double>> computeDFT(const std::vector<double>& signal) {
    int N = signal.size();
    std::vector<std::complex<double>> dft(N);

    for (int k = 0; k < N; ++k) {
        dft[k] = 0;
        for (int n = 0; n < N; ++n) {
            double angle = -2 * PI * k * n / N;
            dft[k] += signal[n] * std::complex<double>(cos(angle), sin(angle));
        }
    }
    return dft;
}

// Алгоритм БПФ (рекурсивная реализация Cooley-Tukey)
std::vector<std::complex<double>> computeFFT(const std::vector<double>& signal) {
    int N = signal.size();

    // Базовый случай
    if (N == 1) {
        return {std::complex<double>(signal[0], 0)};
    }

    // Разделение на четные и нечетные
    std::vector<double> even(N/2), odd(N/2);
    for (int i = 0; i < N/2; ++i) {
        even[i] = signal[2*i];
        odd[i] = signal[2*i + 1];
    }

    // Рекурсивные вызовы
    auto evenFFT = computeFFT(even);
    auto oddFFT = computeFFT(odd);

    // Объединение результатов
    std::vector<std::complex<double>> fft(N);
    for (int k = 0; k < N/2; ++k) {
        double angle = -2 * PI * k / N;
        std::complex<double> twiddle = std::complex<double>(cos(angle), sin(angle));

        fft[k] = evenFFT[k] + twiddle * oddFFT[k];
        fft[k + N/2] = evenFFT[k] - twiddle * oddFFT[k];
    }
    return fft;
}

// Обратное преобразование Фурье
std::vector<double> computeInverseFFT(const std::vector<std::complex<double>>& fft) {
    int N = fft.size();
    std::vector<double> signal(N);

    for (int n = 0; n < N; ++n) {
        std::complex<double> sum = 0;
        for (int k = 0; k < N; ++k) {
            double angle = 2 * PI * k * n / N;
            sum += fft[k] * std::complex<double>(cos(angle), sin(angle));
        }
        signal[n] = sum.real() / N; // Масштабирование
    }
    return signal;
}

// Функция для вычисления погрешности между двумя сигналами
double calculateError(const std::vector<double>& original, const std::vector<double>& reconstructed) {
    double sumSqOriginal = 0.0;
    double sumSqDiff = 0.0;

    for (size_t i = 0; i < original.size(); ++i) {
        sumSqOriginal += original[i] * original[i];
        double diff = original[i] - reconstructed[i];
        sumSqDiff += diff * diff;
    }

    return sqrt(sumSqDiff / sumSqOriginal);
}

// Функция для экспорта данных в файл
void exportToFile(const std::vector<double>& signal, const std::string& filename) {
    std::ofstream outFile(filename);
    for (double value : signal) {
        outFile << value << std::endl;
    }
    outFile.close();
}

void exportComplexToFile(const std::vector<std::complex<double>>& signal, const std::string& filename) {
    std::ofstream outFile(filename);
    for (auto value : signal) {
        outFile << value.real() << " " << value.imag() << std::endl;
    }
    outFile.close();
}

int main() {
    // Начальные параметры для проверки
    int numHarmonics = 3;
    std::vector<double> amplitudes = {1.0, 0.5, 0.3};
    double samplingRate = 1000.0;
    double initialPhase = 0.0;
    int numPoints = 128; // Для БПФ лучше использовать степень двойки
    double deltaT = 1.0 / samplingRate;
    int quantizationLevels = 16;

    // Генерация сигнала
    std::vector<double> analogSignal = generateSignal(numHarmonics, amplitudes, samplingRate, initialPhase, numPoints, deltaT);

    // Преобразование в цифровой сигнал
    std::vector<double> digitalSignal = convertToDigital(analogSignal, quantizationLevels);

    // 1. Прямое преобразование Фурье (DFT)
    auto dftResult = computeDFT(digitalSignal);

    // 2. Быстрое преобразование Фурье (FFT)
    auto fftResult = computeFFT(digitalSignal);

    // 3. Обратное преобразование Фурье для восстановления сигнала
    auto reconstructedSignal = computeInverseFFT(fftResult);

    // 4. Сравнение с оригинальным сигналом
    double error = calculateError(analogSignal, reconstructedSignal);

    std::cout << "Погрешность восстановления: " << error * 100 << "%" << std::endl;

    if (error > ERROR_THRESHOLD) {
        std::cout << "Предупреждение: Погрешность превышает 5%" << std::endl;
    } else {
        std::cout << "Погрешность в допустимых пределах (<= 5%)" << std::endl;
    }

    // Экспорт данных в файлы
    exportToFile(analogSignal, "analog_signal.txt");
    exportToFile(digitalSignal, "digital_signal.txt");
    exportComplexToFile(dftResult, "dft_result.txt");
    exportComplexToFile(fftResult, "fft_result.txt");
    exportToFile(reconstructedSignal, "reconstructed_signal.txt");

    return 0;
}
