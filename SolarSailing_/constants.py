AU = 149597870.7  # км
R_sun = 700000 # км
MU_sun = 132712440019.  # км ^ 3 / с ^2
# Единицы измерения
DU = AU  # Задано в км
TU = 1.  # Задано в секундах; сейчас это 1 день
VU = DU / TU

# Константы
MU_in_units = MU_sun * (TU ** 2 / DU ** 3)  # DU^3/TU^2
AU_in_units = AU / DU
R_sun_in_units = R_sun / DU
