# FFS — pomysły na później

## Detekcja linii: filtr Hessego + Hough z kierunkiem

Zapisane 2026-10-06. Do zrobienia po poprawieniu podstawowego Hougha w `ffs.py`
(siatka ρ/θ, wybór maksimów, maska wejściowa).

### Dlaczego

`line_filter` (jądro pierścieniowe) nie nadaje się jako wejście Hougha: korzysta tylko
z 2 punktów, w których linia przecina pierścień, więc gubi słabe ślady (na klatce
testowej 2% pikseli śladu 1,25σ/px zamiast 19% przy zwykłym gaussie + 3σ) i nie
odrzuca jasnych gwiazd. Dobre jądro do linii powinno sumować sygnał wzdłuż linii, mieć
w poprzek profil ~PSF, sumę zero i — przede wszystkim — odróżniać linię od gwiazdy
przez porównanie odpowiedzi w różnych kierunkach.

### Pomysł A — macierz Hessego w skali śladu („ridge filter”)

- Wygładzić obraz gaussem o σ ≈ szerokość śladu (np. z `frame_fwhm / 2.355`) i policzyć
  drugie pochodne Ixx, Iyy, Ixy (3 rozdzielne sploty, szybkie).
- W każdym pikselu wartości własne λ₁ ≤ λ₂ macierzy [[Ixx, Ixy], [Ixy, Iyy]]:

  | obiekt | λ₁ (w poprzek) | λ₂ (wzdłuż) |
  |---|---|---|
  | jasna linia | silnie ujemna | ≈ 0 |
  | gwiazda | silnie ujemna | silnie ujemna |
  | szum | losowa | losowa |

- Miara liniowości, np. A = (|λ₁| − |λ₂|) / (|λ₁| + |λ₂|) — ~1 linia, ~0 plamka;
  siła = |λ₁|. Wektor własny daje lokalny kierunek linii.
- Uwagi: zbocza gwiazd dają słaby, pierścieniowy sygnał „liniowości”; złe kolumny
  i przelewy też są grzbietami — rozpoznawać osobno (kierunek, nasycenie).
- Literatura: Frangi i in. 1998 (vesselness), Steger 1998 (curvilinear structures).

### Pomysł B — Hough z ograniczeniem kierunku

- Każdy piksel głosuje z wagą z pomysłu A (liniowość × siła), i **tylko na linie
  bliskie swojemu lokalnemu kierunkowi** (±kilka stopni), zamiast na wszystkie 180.
- Efekt: gwiazdy przestają rysować sinusoidy w akumulatorze, szum głosuje dużo
  słabiej, mniej głosów = szybciej. Hough nadal sumuje sygnał wzdłuż całej linii
  („nieskończenie długi filtr dopasowany”), więc słabe ślady zostają.

### Jak ocenić

W `examples/ffs_satellites.ipynb` / na klatkach z `tests/ffs_synth.py`: porównać mapę
wag z pomysłu A z obecną `self.maska` — jaka część pikseli śladu (także słabego, ~1σ/px)
przechodzi, jaka część pikseli gwiazd odpada; potem liczba fałszywych linii i czas
z pomysłem B. Na końcu prawdziwa klatka jk15c (bez satelitów — ma dać 0).

### Dalej (opcjonalnie)

Nir, Zackay & Ofek 2018, „Optimal and Efficient Streak Detection in Astronomical Images”
(AJ) — szybka transformata Radona z poprawną statystyką detekcji smug.

## Przegląd FFS na tle innych bibliotek i lista braków

Zapisane 2026-10-07 (rozmowa przeglądowa, bez zmian w kodzie).

### Pozycja FFS

Żadna pojedyncza funkcja nie przebija najlepszego odpowiednika (photutils/SEP do
detekcji, `Background2D` do tła, GalSim HSM do momentów). Wartość FFS to szybkie QC
klatki w jednym wywołaniu dla obserwatorium + maski dla photutils.

Najmocniejsze/oryginalne elementy:
- `satellites.py`: binomialna istotność względem długości cięciwy (`p·L`), akceptacja
  po medianie segmentów (odporna na gwiazdy), empiryczny szum z równoległych linii,
  klasyfikacja kolumn/wierszy.
- `theta_spread` (statystyka kołowa, okres π) — rzadko spotykana, oddziela prowadzenie/
  wiatr od szumu; zestaw `cpe`, `shape`, `ci` jako szybkie wskaźniki ostrości.

### Kompatybilność z photutils — tarcia

- `x`, `y` to piksele całkowite maksimum, nie centroidy subpikselowe; photutils
  oczekuje `xcentroid`/`ycentroid`.
- Brak obsługi NaN / maski na wejściu (`mk_stats` z NaN daje złe kwantyle).
- Brak wejścia `CCDData`/`NDData`.
- Tło to skalar (mediana), nie mapa tła + mapa rms.
- `saturation = 50000` na sztywno; `box_mag` z umownym punktem zerowym 25.

### Pomysły (od najważniejszych)

1. **Mapa PSF w polu** — FWHM/eliptyczność vs położenie (siatka 3×3 lub wielomian jak
   w `sky_gradient`): tilt, kolimacja, krzywizna pola. Dane są już w `self.stars`.
2. **Centroidy subpikselowe** w tabeli gwiazd (`adaptive_moments` je liczy, nie zapisuje).
3. **Mapa tła + rms** (patrz plan globalnego binowania/tła); potem lokalny próg
   w `find_stars`.
4. Maska nasyconych pikseli i przelewów (kolumny z `satellites.py` → `self.masks`).
5. Maska promieni kosmicznych (astroscrappy lub laplasjan — `laplace_kernel` istnieje).
6. Metoda składająca wszystkie maski w jedną (`|` po `self.masks`).
7. Wejście z NaN/maską; `saturation` jako parametr konstruktora.
8. Niepewności: SNR gwiazdy, błąd strumienia.
9. Zapis wyników do nagłówka FITS (`FWHM`, `ELLIP`, `NSTARS`), eksport regionów DS9.
10. `mk_stats`: `np.partition` lub podpróbka zamiast pełnego sortowania.
11. Pętle po gwiazdach w Pythonie — na razie bez problemu.

### Do sprawdzenia

`mk_stats`: `noise = sqrt(median/gain + rn_noise**2)` jest w ADU tylko gdy `rn_noise`
jest w ADU. Jeśli szum odczytu podawany jest w e⁻, powinno być `(rn_noise/gain)**2`.
Sprawdzić, co przekazuje TOI.

### Przed publikacją

- Usunąć `fwhm_old`, `fwhm_1d_old`, zakomentowane warianty w `cpe`.
- Ujednolicić komentarze/docstringi (angielski).
- Brakujące `@staticmethod` przy `fwhm`, `fwhm_1d`.
- Pozycjonowanie: PyPI (samodzielnie lub w pyaraucaria) + notka RNAAS o detekcji
  satelitów; JOSS dopiero po warstwie masek i tła.
