# Fortran CSV Reader


## Verzeichnisstruktur

```
.
├── build.sh                # Skript zum Kompilieren der Fortran-Bibliothek
├── fpm.toml                # Fortran Package Manager Konfiguration
├── python/
│   └── read_csv.py         # Python-Wrapper und Testskript
├── R/
│   └── read_csv.R          # R-Wrapper und Testskript
├── src/
│   └── csv_module/         # Fortran-Quellcode-Module
│       ├── csv_parser_module.f90
│       ├── csv_file_reader_module.f90
│       └── csv_read_*.f90
└── test/
    ├── asserts.f90         # Test-Hilfsfunktionen
    ├── run_tests.f90       # FPM-Test-Runner
    └── mod_test_*.f90      # Test-Suiten für jedes Modul
```

---

## Build-Anweisungen

Die gesamte Fortran-Bibliothek wird mit dem `build.sh`-Skript kompiliert. Dieses Skript verwendet `fpm`, um die Objektdateien zu erstellen, und linkt sie dann manuell zu einer Shared Library.

**Voraussetzungen:** Ein Fortran-Compiler (z.B. `gfortran`) und `fpm`.

**Befehl:**
Führen Sie im Stammverzeichnis des Projekts den folgenden Befehl aus:
```bash
./build.sh
```
Nach erfolgreicher Ausführung befindet sich die kompilierte Shared Library unter `build/libtensor-omics.so`. Dieser Schritt ist **erforderlich**, bevor der Python- oder R-Wrapper verwendet werden kann.

---

## Verwendung

Die Bibliothek ist so konzipiert, dass die aufrufende Sprache (Python, R, etc.) den Leseprozess steuert. Der typische Arbeitsablauf ist immer:

1.  **Dimensionen abfragen:** Rufen Sie `get_csv_dims_*` auf, um die Anzahl der Zeilen und Spalten zu ermitteln.
2.  **Speicher allozieren:** Erstellen Sie in Ihrer Sprache (z.B. mit NumPy) leere Arrays mit den korrekten Dimensionen, um die Ergebnisse zu speichern.
3.  **Daten als Text lesen:** Rufen Sie `read_csv_to_strings_*` auf, um die gesamte CSV-Datei in ein 2D-String-Array (bzw. ein Array von ASCII-Codes) zu lesen.
4.  **Spalten konvertieren:** Rufen Sie die spezifischen `read_*_columns_*`-Funktionen für jede gewünschte Spalte auf, um die Textdaten in die Zieldatentypen zu konvertieren.

### Python-Verwendung

Das Skript `python/read_csv.py` enthält alle notwendigen Funktionen und ein vollständiges Testbeispiel.

**Beispiel:**
```python
# Annahme: Sie befinden sich im python/ Verzeichnis
from read_csv import setup_fortran_functions, _load_fortran_library
import numpy as np
import ctypes

# 1. Bibliothek laden und Funktionen einrichten
fortran_lib = _load_fortran_library()
funcs = setup_fortran_functions(fortran_lib)

# 2. Dateinamen vorbereiten
filename = "../data.csv" # Beispielpfad
fname_ascii = np.array([ord(c) for c in filename], dtype=np.int32)

# 3. Dimensionen abfragen
num_rows, num_cols, status = ctypes.c_int(), ctypes.c_int(), ctypes.c_int()
funcs['get_csv_dims'](fname_ascii, len(filename), True, ord(','),
                      ctypes.byref(num_rows), ctypes.byref(num_cols), ctypes.byref(status))

nr, nc = num_rows.value, num_cols.value

# 4. Speicher für String-Daten allozieren und Daten lesen
# (Dieser Teil wird oft von einer Hilfsfunktion gekapselt)
# ... siehe read_csv_to_strings in wrapper.py für das vollständige Beispiel

# 5. Gewünschte Spalten konvertieren
# ... siehe read_integer_columns etc. in wrapper.py für das vollständige Beispiel
```

### R-Verwendung

Das Skript `R/read_csv.R` dient als vollständiges Beispiel.

**Beispiel:**
```R
# 1. Bibliothek laden
lib_path <- file.path("..", "build", "libtensor-omics.so")
dyn.load(lib_path)

# 2. Dateinamen vorbereiten
filename <- "../data.csv" # Beispielpfad
fname_ascii <- as.integer(charToRaw(filename))

# 3. Dimensionen abfragen
res_dims <- .Fortran("get_csv_dims_r",
                     filename_ascii = fname_ascii,
                     fn_len = as.integer(length(fname_ascii)),
                     has_header = TRUE,
                     delimiter_ascii = as.integer(charToRaw(",")),
                     num_rows = integer(1),
                     num_cols = integer(1),
                     status = integer(1))

nr <- res_dims$num_rows
nc <- res_dims$num_cols

# 4. Speicher für String-Daten allozieren und Daten lesen
# ... siehe test_csv_reader.R für das vollständige Beispiel

# 5. Gewünschte Spalten konvertieren
# ... siehe test_csv_reader.R für das vollständige Beispiel
```

### Fortran-Verwendung (Beispiel)

Ein reines Fortran-Programm kann die Module direkt verwenden und linken.

```fortran
PROGRAM example
    USE csv_file_reader_module
    USE csv_read_int_moduleays für die Roh-String-Daten und für jede typisierte Sp
    ! ... weitere USE-Anweisungen ...
    IMPLICIT NONE
    
    CHARACTER(LEN=64) :: filename
    CHARACTER(LEN=512), ALLOCATABLE :: data_str(:,:)
    INTEGER(INT32), ALLOCATABLE :: int_col(:,:)
    ! ... weitere Deklarationen ...
    
    filename = 'my_data.csv'
    
    ! Lese CSV als Strings
    CALL read_csv_to_strings(filename, .TRUE., ',', header, data_str, status)
    
    ! Konvertiere die erste Spalte zu Integer
    CALL read_integer_columns(data_str, [1], int_col, status)
    
    ! Gib die Integer-Spalte aus
    PRINT *, int_col
    
END PROGRAM example
```
**Kompilierungsbefehl:**
```bash
# Finde zuerst das .mod-Verzeichnis
MOD_DIR=$(find build -name "*.mod" -printf "%h" | head -n 1)
# Kompiliere und linke
gfortran example.f90 -I"$MOD_DIR" -Lbuild -ltensor-omics -o example
```

---

## API-Referenz

| Funktionalität | Fortran-Subroutine | C-Wrapper | R-Wrapper |
| :--- | :--- | :--- | :--- |
| CSV-Dimensionen ermitteln | (intern) | `get_csv_dims_c` | `get_csv_dims_r` |
| CSV als Text einlesen | `read_csv_to_strings` | `read_csv_to_strings_c` | `read_csv_to_strings_r` |
| Zeile parsen | `parse_line` | `parse_line_c` | `parse_line_r` |
| Integer-Spalten lesen | `read_integer_columns` | `read_integer_columns_c` | `read_integer_columns_r` |
| Real-Spalten lesen | `read_real_columns` | `read_real_columns_c` | `read_real_columns_r` |
| Logical-Spalten lesen | `read_logical_columns`| `read_logical_columns_c`| `read_logical_columns_r`|
| Character-Spalten lesen | `read_character_columns`|`read_character_columns_c`|`read_character_columns_r`|
| Complex-Spalten lesen | `read_complex_columns`|`read_complex_columns_c` |`read_complex_columns_r` |

