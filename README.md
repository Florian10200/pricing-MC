# **PRICER, ENSAE 23-24**

Ce code permet de déterminer les prix d'un large évantail d'options (européennes, barrières, lookback et asiatiques).

Il s'appuie sur le modèle de Black-Scholes-Merton, et dans le cas où des formules explicites ne s'appliquent pas, les prix sont calculés par méthode de Monte-Carlo.

Pour exécuter le code, il est nécessaire de disposer de la version de Code::Blocks v20.03 ou d'une version postérieure.

# Structure

La structure de notre code est articulée autour de trois fichiers principaux. Tout d'abord, le fichier pricer.h recense toutes les fonctions et classes indispensables à la réalisation du projet. Ensuite, le fichier pricer.cpp constitue l'essence même de notre implémentation, hébergeant l'ensemble des fonctions de pricing développées, il contient le code source. Enfin, le main.cpp, offre la possibilité de générer le prix de l'option en fonction de ses caractéristiques spécifiques, intégrant ainsi une dimension opérationnelle à notre travail.

Vous pouvez tester le code avec les exemples que nous avons implémentés, ou modifier les paramètres des options pour calculer d'autres prix.

# Auteurs

Florian LAVA

Hugues MATON

Taddéo VIVET
