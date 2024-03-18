# -*- coding: utf-8 -*-
glagolitsa = "А,Б,В,Г,Д,Е,Ё,Ж,З,И,Й,К,Л,М,Н,О,П,Р,С,Т,У,Ф,Х,Ц,Ч,Ш,Щ,Ъ,Ы,Ь,Э,Ю,Я"
Glagolitsa = glagolitsa.split(',')
names = ["Денисов Артемий",
         "Ильин Вячеслав",
         "Крупский Андрей",
         "Тарасевич Артем",
         "Шибанов Игорь",
         "Архангельский Андрей",
         "Багринцев Михаил",
         "Загребин Григорий",
         "Петрова Елизавета",
         "Агафонов Вадим",
         "Агафонов Дадим",
         "Шинкарева Алена"]


def get_index(ch):
    return Glagolitsa.index(ch)


ord_A = get_index('А')

for name in names:
    lastname, firstname = name.split()
    var = (get_index(lastname[0]) - ord_A + get_index(firstname[0]) - ord_A) % 2 + 1
    print(var, ": ", name)


for letter in Glagolitsa:
    print(ord(letter), end=" ")
