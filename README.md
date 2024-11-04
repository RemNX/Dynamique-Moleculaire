# Atomistique

Bienvenue dans le projet **Atomistique** ! Ce dépôt contient des implémentations et des démonstrations liées à l'atomistique, notamment des méthodes de simulation et d'analyse des systèmes atomiques.

## Table des matières

- [Présentation](#présentation)
- [Fonctionnalités](#fonctionnalités)
- [Utilisation](#utilisation)
- [Contact](#contact)

## Présentation

Le projet Atomistique a pour objectif d'explorer les concepts fondamentaux de la dynamique moléculaire et de fournir des outils pour simuler des systèmes atomiques. Nous avons développé ce code dans le cadre du cours Simulation atomistique des matériaux à la faculté des sciences de Montpellier. Nous travaillons sur un système de N particules en interaction Lennard-Jones dans une boîte en 2D, avec des conditions aux bords périodiques. Nous avons deux codes principaux en C : l'un qui vérifie l'exactitude de l'implémentation des algorithmes, l'autre qui calcule les grandeurs thermodynamiques à l'aide des relations de la physique statistique. Un code en Python, `anime.py`, récupère les données concernant les positions des particules afin de réaliser une animation en utilisant FuncAnimation de Matplotlib.

## Fonctionnalités

- Simulation de systèmes Lennard-Jones.
- Implémentations de différentes méthodes numériques, y compris la méthode d'Euler et la méthode de Velocity-Verlet.
- Vérification d'algorithmes.
- Animation visuelle des particules dans une boîte.
- Analyse et visualisation des résultats.
- Fusion de cristaux et inclusion de défauts.

## Utilisation 

- Prenez le code `main_simulation.c` et choisissez les paramètres que vous voulez, comme le nombre de particules et la taille de la boîte.
- Lancez le code avec `cc main_simulation.c -lm && ./a.out && python3 anime.py` si vous êtes sur terminal.
- Pour plus d'informations, suivez les commentaires présents dans le fichier `main_simulation.c`.

## Contact

Si vous avez des questions, des problèmes ou des suggestions concernant le projet, n'hésitez pas à ouvrir une issue sur ce dépôt pour signaler un problème ou poser une question.

