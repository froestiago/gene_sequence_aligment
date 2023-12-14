
# algoritmo Needleman–Wunsch
class NW(object):
    _secuencia1: [str]
    _secuencia2: [str]
    _matriz = []
    _match: int
    _mismatch: int
    _gap: int
    _costo_ext: int

    _resultados = {}

    def __init__(self, seq1, seq2, match, mismatch, gap, costo_ext):
        """ Constructor de la clase

        Args:
            seq1 (str): Secuencia 1
            seq2 (str): secuencia 2
            match (int): Valor coincidencia
            mismatch (int): Costo falta de coincidencia
            gap (int): Costo gap
            costo_ext (int): Costo extension

        """

        self._secuencia1 = ["-"] + list(seq1) # m +1
        self._secuencia2 = ["-"] + list(seq2) # n +1
        self._match = match
        self._mismatch = mismatch
        self._gap = gap
        self._costo_ext = costo_ext

        self._rellenar_matriz()
        self._backtracking()
        self._get_identidad_similitud()
        self._imprimir_resultados()
        # self.imprimir_matriz()

    def _rellenar_matriz(self):
        """ Rellena la matriz mientra calcula los valores """
        for i in range(len(self._secuencia1)):
            print(f'{i}/{len(self._secuencia1)}')
            self._matriz.append([])
            for j in range(len(self._secuencia2)):
                self._matriz[i].append(self._calcular_valor(i, j))

    def _calcular_valor(self, i, j):
        """ Calcula el valor de la posicion i, j de la matriz 
        
        Args:
            i (int): Posicion fila
            j (int): Posicion columna

        Returns:
            int: Valor de la posicion i, j
        """

        if i == 0 and j == 0:
            return 0
        elif i == 0:
            return self._matriz[i][j-1] + self._gap
        elif j == 0:
            return self._matriz[i-1][j] + self._gap

        else:
            diagonal = None
            if self._secuencia1[i] == self._secuencia2[j]:
                diagonal = self._matriz[i-1][j-1] + self._match
            else:
                diagonal = self._matriz[i-1][j-1] + self._mismatch

            return max(diagonal, self._matriz[i-1][j] + self._gap, self._matriz[i][j-1] + self._gap)

    def _backtracking(self):
        """ Realiza el backtracking de la matriz """

        seq1 = str()
        seq2 = str()
        extension = 0
        i = len(self._matriz) - 1
        j = len(self._matriz[0]) - 1
        while i > 0 or j > 0:
            if i == 0:
                seq1 = "-" + seq1
                seq2 = self._secuencia2[j] + seq2
                j -= 1
                extension += 1

            elif j == 0:
                seq1 = self._secuencia1[i] + seq1
                seq2 = "-" + seq2
                i -= 1
                extension += 1

            elif self._secuencia1[i] == self._secuencia2[j]:
                # diagonal
                seq1 = self._secuencia1[i] + seq1
                seq2 = self._secuencia2[j] + seq2
                i -= 1
                j -= 1

            elif self._matriz[i-1][j-1] >= max(self._matriz[i-1][j], self._matriz[i][j-1]):
                # diagonal                      arriba                izquierda
                seq1 = self._secuencia1[i] + seq1
                seq2 = self._secuencia2[j] + seq2
                i -= 1
                j -= 1

            elif self._matriz[i-1][j] >= max(self._matriz[i-1][j-1], self._matriz[i][j-1]):
                # arriba                        diagonal             izquierda
                seq1 = self._secuencia1[i] + seq1
                seq2 = "-" + seq2
                i -= 1
                extension += 1
            elif self._matriz[i][j-1] >= max(self._matriz[i-1][j-1], self._matriz[i-1][j]):
                # izquierda                     diagonal             arriba
                seq1 = "-" + seq1
                seq2 = self._secuencia2[j] + seq2
                j -= 1
                extension += 1

        self._resultados["secuencia1"] = seq1
        self._resultados["secuencia2"] = seq2
        self._resultados["costo extension"] = self._costo_ext * extension

    def _get_identidad_similitud(self):
        """ Calcula estadísticas """
        identidad = 0
        similitud = 0
        gaps = 0
        largo = len(self._resultados["secuencia1"])

        for i in range(len(self._resultados["secuencia1"])):
            if self._resultados["secuencia1"][i] == self._resultados["secuencia2"][i]:
                identidad += 1
            elif self._resultados["secuencia1"][i] != "-" and self._resultados["secuencia2"][i] != "-":
                similitud += 1
            elif self._resultados["secuencia1"][i] == "-" or self._resultados["secuencia2"][i] == "-":
                gaps += 1

        self._resultados["coincidencias"] = identidad
        self._resultados["similitudes"] = similitud
        self._resultados["gaps"] = gaps
        self._resultados["score"] = (identidad * self._match) + (similitud * self._mismatch) + (self._costo_ext * gaps)
        self._resultados["costo_extension_total"] = self._costo_ext * gaps
        self._resultados["porcentaje_identidad"] = f"{(identidad / largo)*100}%"
        self._resultados["porcentaje_similitud"] = f"{((similitud) / largo)*100}%"
        self._resultados["porcentaje_gaps"] = f"{(gaps / largo)*100}%"
        

    def _imprimir_resultados(self):
        """ Imprime los resultados """
        for key, value in self._resultados.items():
            print(f"{key} : {value}")

    def imprimir_matriz(self):
        """ Imprime la matriz """
        print(f"{'':>8}{'    '.join(self._secuencia2):>4}")
        for i in range(len(self._matriz)):
            print(f"{self._secuencia1[i]:>4}", end=" ")
            for j in range(len(self._matriz[i])):
                print(f"{self._matriz[i][j]:>4}", end=" ")
            print()
