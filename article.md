---
title : "New_simulator"
date : "2022-06-27T15:52:47+02:00"
author : "Daniel Simu & Ilse Arwert"
authorTwitter : "" #do not include @
cover : ""
tags : []
keywords : []
description : ""
showFullContent : false
readingTime : true
draft : true
toc : true
---

# Introduction

Siteswap simulators and siteswap pattern generators are very useful tools to jugglers. Countless tricks have been discovered with the help of generators, and tricks have been learnt with the help of simulators or animators.

Once you understand a bit of siteswap it is a lot of fun to plug in various sequences of numbers into a simulator and test all kinds of patterns to see what they would look like. But, there are a lot of features of juggling that can't be described with siteswaps, and almost all juggling simulators that I'm aware of are siteswap based.

[Rhythmic Catches](../rhythmic_catches/) and [IMBO](../improved_body_trick_notation/) have been designed with trick generation and simulation mind. Together with Ilse Arwert I attempted to modify Juggling Lab to make it suitable for these notation systems.

# Why?

What are the possible motivations to generate or simulate tricks?

## Generation

Trick generation allows you to explore trick discovery in a different way than by testing out ideas with props in a studio. It can be faster and easier, and with a good tool you can filter for exactly the kind of trick you are looking for. For example you might look for tricks that take a certain amount of beats, involve a certain amount of jugglers, include certain kinds of throws or exclude throws above a certain height. Also one can look for variations on a particular trick. With a more detailed notation system, more parameters become available to discover or to filter for.

There is also a more philosophical aspect to trick generation. If a trick is generated by a computer, rather than organically found by a juggler, who is then the inventor of that trick? Can one still claim ownership of an idea if that idea is the natural result of a model of juggling that is able to generate that idea? I hope that with better generators we can further this debate.

## Simulation

Once you have a notation for a trick, you can try and visualize it, or simply learn it and perform it, to see what the notation means. However, a juggling animator or simulator is able to do so much more swiftly. From the animation you can more easily understand what goes on, decide if you like it or wish to modify it. You can slow it down and pause it, which can help a lot with comprehending the pattern.

There exist also juggling simulators that allow you to edit patterns, so that once you have the simulation you can change the timings, the prop types, or the positions of hands and bodies. Through this editing you can manually create completely new patterns.

# Our own simulator

Inspired by these various motivations, and excited to see the true power of Rhythmic Catches and IMBO in action, Ilse Arwert and myself set out to program a juggling simulator that could handle some of the features of Rhythmic Catches. 

We chose to do this by modifying the open source software [Juggling Lab](https://jugglinglab.org/), for which we had full support and encouragement of its creator Jack Boyce.

## Why modify Juggling Lab?

When we considered to build our own simulator, we had the option to either start from scratch or build upon existing software. To save time we chose to edit existing software, and for this we considered two options.

Juggling lab is open source software written in Java. It has been developed [since at least 1997](http://juggling.org/programs/java/JuggleAnim/ja.html) (then called JuggleAnim) and is still being [updated regularly](https://github.com/jkboyce/jugglinglab). It has many different features, for example the ability to control hand positions which is useful if we wanted to simulate body throws.

The other option was JoePass!, which is open source software written in C++. However, it features less intuitive options to move the limbs of the jugglers, and also I was simply less familiar with it. For this reason, we chose to work with Juggling Lab.

As a bonus, I was already in touch with Jack Boyce, who encouraged us to make changes to the software.

## Original goals 

It would have been nice if we could simulate all of Rhythmic Catches and IMBO, but this was a larger task than we imagined to be possible for ourselves within the limited time we wanted to spend on this project.

Therefor we set out the following targets:

- Adding legs to the juggler in Juggling Lab
- Allowing balls to be tossed and caught from any position on the body
- Extending the internal xml language called JML to be able to access these new features.

This would allow me to later write a new interface that could translate Rhythmic Catches and IMBO into extended JML and make it do what I wanted it to do.

## JML is it good? v4 challenges

## Challenges with legacy code, specifically jugglinglab

## What we achieved, Screenshots & link to github repo

# recommendations of how one should build the simulator of the future

# Appendix: List of juggling software

Here is a short list of interesting juggling software I came across whilst researching this project:

- [Juggling Lab](https://jugglinglab.org/), 3D multi feature juggling simulator, pattern generator, pattern editor, especially useful for ball juggling.
- [JoePass!](http://westerboer.net/w/?page_id=151), 3D multi feature juggling simulator, especially useful for passing.
- [Passist](https://dev.passist.org/), passing siteswap generator, diagram generator and animator
- [Prechacthis.org](http://www.prechacthis.org/), prechac pattern generator and editor
- [Juggle Suggest 2](https://joshmermelstein.com/juggle-suggest2/), siteswap autocompletion and animation
- [Polyrhythmic Fountain](https://joshmermelstein.com/polyrhythmic-fountain/), polyrhythmic juggling simulator
- [Juggling Toolbox](https://www.jonglage.net/jugglingTB.html), Generators for various different notation systems
- [Power Juggler](http://vcg.isti.cnr.it/~tarini/index.php?3?juggling.html), simulator for negative siteswaps/antimatter.
- Realistic Juggling Simulator by Pedro Teodoro, juggling siteswaps with imperfect accuracy. I can't find an active link, but a copy can be found in the Juggling Toolbox.

Of course this list is far far far from extensive, it's just some stuff that I find particularly interesting.