#include <vector>
#include <iomanip>
#include <iostream>
#include <random>
#include <omp.h>
#include <chrono>
#include <functional>

bool isSorted(const std::vector<int>& arr) {
    for (size_t i = 1; i < arr.size(); ++i) {
        if (arr[i - 1] > arr[i]) {
            std::cout << "ERROR: ARRAY IS NOT SORTED";
            return false; // Найдено нарушение порядка, массив не отсортирован
        }
    }
    return true; // Массив отсортирован
}

#pragma region merge sort
// Последовательная сортировка слиянием
void merge(std::vector<int>& arr, int left, int mid, int right) {
    int n1 = mid - left + 1;
    int n2 = right - mid;
    std::vector<int> L(n1), R(n2);

    for (int i = 0; i < n1; ++i)
        L[i] = arr[left + i];
    for (int i = 0; i < n2; ++i)
        R[i] = arr[mid + 1 + i];

    int i = 0, j = 0, k = left;
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            ++i;
        }
        else {
            arr[k] = R[j];
            ++j;
        }
        ++k;
    }

    while (i < n1) {
        arr[k] = L[i];
        ++i;
        ++k;
    }

    while (j < n2) {
        arr[k] = R[j];
        ++j;
        ++k;
    }
}

void mergeSort(std::vector<int>& arr, int left, int right) {
    if (left < right) {
        int mid = left + (right - left) / 2;
        mergeSort(arr, left, mid);
        mergeSort(arr, mid + 1, right);
        merge(arr, left, mid, right);
    }
}

// Параллельная сортировка слиянием
void parallelMergeSort(std::vector<int>& arr, int left, int right, int depth) {
    if (left < right) {
        if (depth <= 0) {
            mergeSort(arr, left, right);
        }
        else {
            int mid = left + (right - left) / 2;

            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    parallelMergeSort(arr, left, mid, depth - 1);
                }
                #pragma omp section
                {
                    parallelMergeSort(arr, mid + 1, right, depth - 1);
                }
            }
            merge(arr, left, mid, right);
        }
    }
}
#pragma endregion


#pragma region buble sort
// Функция последовательной сортировки пузырьком
void bubbleSort(std::vector<int>& arr) {
    int n = arr.size();
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - i - 1; ++j) {
            if (arr[j] > arr[j + 1]) {
                std::swap(arr[j], arr[j + 1]);
            }
        }
    }
    isSorted(arr);
}

// Функция паралельным сортировки пузырьком
void parallelBubbleSort(std::vector<int>& arr) {
    int n = arr.size();
    bool sorted = false;

    while (!sorted) {
        sorted = true;

        // Четные фазы
        #pragma omp parallel for shared(arr, sorted)
        for (int i = 1; i < n - 1; i += 2) {
            if (arr[i] > arr[i + 1]) {
                std::swap(arr[i], arr[i + 1]);
                sorted = false;
            }
        }

        // Нечетные фазы
        #pragma omp parallel for shared(arr, sorted)
        for (int i = 0; i < n - 1; i += 2) {
            if (arr[i] > arr[i + 1]) {
                std::swap(arr[i], arr[i + 1]);
                sorted = false;
            }
        }
    }
    isSorted(arr);
}
#pragma endregion

#pragma region quick sort
// Функция для разделения массива
int partition(std::vector<int>& arr, int low, int high) {
    int pivot = arr[high];
    int i = low - 1;

    for (int j = low; j < high; ++j) {
        if (arr[j] < pivot) {
            ++i;
            std::swap(arr[i], arr[j]);
        }
    }
    std::swap(arr[i + 1], arr[high]);
    return i + 1;
}

// Последовательная быстрая сортировка
void quickSort(std::vector<int>& arr, int low, int high) {
    if (low < high) {
        int pi = partition(arr, low, high);
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}

// Параллельная быстрая сортировка
void parallelQuickSort(std::vector<int>& arr, int low, int high, int depth) {
    if (low < high) {
        int pi = partition(arr, low, high);

        if (depth > 0) {
            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    parallelQuickSort(arr, low, pi - 1, depth - 1);
                }
                #pragma omp section
                {
                    parallelQuickSort(arr, pi + 1, high, depth - 1);
                }
            }
        }
        else {
            parallelBubbleSort(arr);
        }
    }
}
#pragma endregion

std::vector<int> generateArray(int min, int max, int N) {
    std::vector<int> array(N);
    std::random_device rd;  // Источник случайности
    std::mt19937 gen(rd()); // Генератор случайных чисел
    std::uniform_int_distribution<> dis(min, max); // Распределение от min до max

    for (int i = 0; i < N; ++i) {
        array[i] = dis(gen);
    }

    return array;
}

// Функция для сравнения и вывода результатов сортировки слиянием
void compareMergeSort(std::vector<int>& arraySizes, std::vector<int>& threadCounts) {
    std::cout << "\n\nSize\tSequential Merge\tParallel Merge (2)\tParallel Merge (4)\tParallel Merge (8)\n";

    for (int size : arraySizes) {
        // Генерация массива
        std::vector<int> arr = generateArray(1, 1000, size);

        // Последовательная сортировка слиянием
        std::vector<int> arrSeq = arr; // Копия массива
        auto start = std::chrono::high_resolution_clock::now();
        mergeSort(arrSeq, 0, arrSeq.size() - 1);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> durationSeqMerge = end - start;

        // Проверка, отсортирован ли массив
        isSorted(arrSeq);

        std::cout << std::setw(8) << size << "\t" << std::fixed << std::setprecision(6) << durationSeqMerge.count() << "\t\t";

        // Параллельная сортировка слиянием
        for (int threads : threadCounts) {
            std::vector<int> arrPar = arr; // Копия массива
            start = std::chrono::high_resolution_clock::now();
            omp_set_num_threads(threads);
            parallelMergeSort(arrPar, 0, arrPar.size() - 1, threads);
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> durationParMerge = end - start;

            // Проверка, отсортирован ли массив
            isSorted(arrPar);

            std::cout << std::setw(8) << std::fixed << std::setprecision(6) << durationParMerge.count() << "\t\t";
        }
        std::cout << "\n";
    }
}

void compareBubbleSort(std::vector<int>& arraySizes, std::vector<int>& threadCounts) {
    std::cout << "Size\tSequential Bubble\tParallel Bubble (2)\tParallel Bubble (4)\tParallel Bubble (8)\n";

    for (int size : arraySizes) {
        // Генерация массива
        std::vector<int> arr = generateArray(1, 1000, size);

        // Последовательная пузырьковая сортировка
        auto start = std::chrono::high_resolution_clock::now();
        bubbleSort(arr);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> durationSeqBubble = end - start;

        std::cout << std::setw(8) << size << "\t" << std::fixed << std::setprecision(6) << durationSeqBubble.count() << "\t\t";
        for (int threads : threadCounts) {
            start = std::chrono::high_resolution_clock::now();
            omp_set_num_threads(threads);
            std::vector<int> arrParBubble = arr; // Копия массива
            parallelBubbleSort(arrParBubble);
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> durationParBubble = end - start;
            std::cout << std::setw(8) << std::fixed << std::setprecision(6) << durationParBubble.count() << "\t\t";
        }
        std::cout << "\n";
    }
}

void compareQuickSort(std::vector<int>& arraySizes, std::vector<int>& threadCounts) {
    std::cout << "\n\nSize\tSequential Quick\tParallel Quick (2)\tParallel Quick (4)\tParallel Quick (8)\n";

    for (int size : arraySizes) {
        // Генерация массива
        std::vector<int> arr = generateArray(1, 1000, size);

        // Последовательная быстрая сортировка
        auto start = std::chrono::high_resolution_clock::now();
        quickSort(arr, 0, arr.size() - 1);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> durationSeqQuick = end - start;

        std::cout << std::setw(8) << size << "\t" << std::fixed << std::setprecision(6) << durationSeqQuick.count() << "\t\t";
        for (int threads : threadCounts) {
            start = std::chrono::high_resolution_clock::now();
            omp_set_num_threads(threads);
            std::vector<int> arrParQuick = arr; // Копия массива
            parallelQuickSort(arrParQuick, 0, arrParQuick.size() - 1, threads);
            isSorted(arrParQuick);
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> durationParQuick = end - start;
            std::cout << std::setw(8) << std::fixed << std::setprecision(6) << durationParQuick.count() << "\t\t";
        }
        std::cout << "\n";
    }
}



int main() {
    std::vector<int> threadCounts = { 2, 4, 8 }; // Количество потоков
    std::vector<int> arraySizes = { 100, 1000, 10000 }; // Размеры массивов

    compareBubbleSort(arraySizes, threadCounts);
    arraySizes = { 100, 1000, 10000, 100000, 1000000 }; // Размеры массивов

    compareQuickSort(arraySizes, threadCounts);

    //threadCounts = { 2, 4, 8 }; // Количество потоков
    //arraySizes = { 100, 1000, 10000, 100000, 1000000 }; // Размеры массивов
    //compareMergeSort(arraySizes, threadCounts);
    //compareQuickSort(arraySizes, threadCounts);

    return 0;
}


