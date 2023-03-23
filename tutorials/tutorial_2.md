Friday 2023-03-24 09:00

Functions

## Basics

What are they?

> A function is a block of code that you can use over and over again. It's like
> a recipe that you can follow to make something. In programming, a function
> takes some input, processes it, and returns some output. Think of it as a
> machine that takes in some raw material, processes it, and produces a
> finished product. Functions are used to make your code more modular and
> organized. Instead of writing the same code over and over again, you can
> write a function to do it for you. This makes your code easier to read,
> maintain, and debug.

How do you write them?

When do you use them?
Why do you use them?

What is an argument?
Can we have a function without an argument?
What if we don't want to put a common argument value in? -> default functions

What is a good docstring for a function?

```python
def area(width, height):
    """Calculate the area of a rectangle

    Arguments:
    width: float -- the width of the rectangle
    height: float -- the height of the rectangle

    Returns:
    float -- the area of the rectangle
    """
    return width * height
```

### Scope

Function scope is really important. Do broccoli to show what is in and outside of the scope.

## A problem: a list of books


Titles:

```
De wetten
Darwinian Populations
The Social Construction of What?
The Right to Sex
Enactive Psychiatry
```

Authors:

```
Connie Palmen
Peter Godfrey-Smith
Ian Hacking
Amia Srinivasan
Sanneke de Haan
```

Years:

```
1991
2009
1999
2021
2021
```

## Printing the list

How can we make a list like this:

```
1. Connie Palmen, De wetten (1991)
2. Peter Godfrey-Smith, Darwinian Populations (2009)
3. ...
```

First, how can we print the number (`1.`) followed by the book title?
What would the function look like that does that for us?
What would a good docstring be?

Can we use a function that gives us the number and the book title as a list of pairs?
What would a good docstring be?

Let's integrate that. How?

How can we list both the author name and the book title?
What would a good docstring be?

How do we need to modify the original function?

How can we also add the year?
What would a good docstring be?

Can we write the year in parentheses?
Is it nice to use a seperate function for that?
When is that the case, and when not?
What would a good docstring be?


Okay, what if we want to cut off the books past a certain age?
What would a good docstring be?
