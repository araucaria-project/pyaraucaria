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
