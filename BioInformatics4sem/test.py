from sympy import series, symbols, Function, Eq, solve, expand

# ty chatgpt

# 1

# Кажется понятно. Ладно

t = symbols('t')

A = (1 - t - t ** 2) / (1 + 2 * t ** 3)

series_expansion = series(A, t, n=11).removeO()

a_10 = series_expansion.coeff(t, 10)
print(a_10)

# 2

# Ручками легче в итоге

# 3

r = symbols('r')

# Пример
# an = 3an-1 + 6an-2 - 8an-3 is
# r^n = 3r^(n-1) + 6r^(n-2) - 8r^(n-3)
# r^3 - 3r^2 - 6r + 8 = 0

# Ладно

char_eq = r ** 3 - 3 * r ** 2 + 6 * r + 8

roots = solve(char_eq, r)

root_types = [(root.evalf(), root.is_real) for root in roots]

print(root_types)

# 4

# Для сета - (1 + t ^ вес) ^ кол-во
# Для мультисета - (1 / (1 - t ^ {вес})) ^ кол-во

# Set+ = set - 1
# Mset+ = Mset - 1

# Seq = 1 / (1 - A(t))
# Seq+ = A(t) / (1 - A(t))

# В примере 2 объекта по 1 и два по два

t = symbols('t')

# Формула по весам для сета
# Для мультсета можно без expand (требуется подтверждение)

G_PSet_X = expand((1 + t) ** 2 * (1 + t ** 2) ** 2) - 1
G_MSet_X = (1 / (1 - t)) ** 2 * (1 / (1 - t ** 2)) ** 2 - 1

G_Seq_PSet_X = G_PSet_X / (1 - G_PSet_X)
print(G_Seq_PSet_X)

# 5

# Ручками