# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 15:46:30 2016

@author: Alonyan
"""
import numpy as np
import scipy.io as sio


def num2str(num, precision):
    return "%0.*f" % (precision, num)


def colorcode(datax, datay):
    from scipy import interpolate
    import numpy as np

    H, xedges, yedges = np.histogram2d(datax, datay, bins=30)
    xedges = (xedges[:-1] + xedges[1:]) / 2
    yedges = (yedges[:-1] + yedges[1:]) / 2
    f = interpolate.RectBivariateSpline(xedges, yedges, H)

    z = np.array([])
    for i in datax.index:
        z = np.append(z, f(datax[i], datay[i]))
    # z=(z-min(z))/(max(z)-min(z))
    z[z < 0] = 0
    idx = z.argsort()
    return z, idx


class kmeans:
    def __init__(self, X, K):
        # Initialize to K random centers
        oldmu = X.sample(K).values  # np.random.sample(X, K)
        mu = X.sample(K).values  # np.random.sample(X, K)
        while not _has_converged(mu, oldmu):
            oldmu = mu
            # Assign all points in X to clusters
            clusters = _cluster_points(X, mu)
            # Reevaluate centers
            mu = _reevaluate_centers(oldmu, clusters)
        self.mu = mu
        self.clusters = clusters
        # return(mu, clusters)


def _cluster_points(X, mu):
    clusters = {}
    for x in X:
        bestmukey = min(
            [(i[0], np.linalg.norm(x - mu[i[0]])) for i in enumerate(mu)],
            key=lambda t: t[1],
        )[0]
        try:
            clusters[bestmukey].append(x)
        except KeyError:
            clusters[bestmukey] = [x]
    return clusters


def _reevaluate_centers(mu, clusters):
    newmu = []
    keys = sorted(clusters.keys())
    for k in keys:
        newmu.append(np.mean(clusters[k], axis=0))
    return newmu


def _has_converged(mu, oldmu):
    return set(mu) == set(oldmu)


def makeTicks():
    a = np.outer(np.arange(1, 10), 10 ** np.arange(1, 2)).T.reshape((1, -1)).squeeze()
    ticks = np.append(-a[::-1], 0)
    ticks = np.append(-100, ticks)
    a = np.outer(np.arange(1, 10), 10 ** np.arange(1, 6)).T.reshape((1, -1)).squeeze()
    ticks = np.append(ticks, a[:])
    emptvec = ["", "", "", "", "", "", "", ""]
    ticklabels = (
        ["-0.1"]
        + emptvec
        + [""]
        + ["0"]
        + emptvec
        + [""]
        + ["0.1"]
        + emptvec
        + ["1"]
        + emptvec
        + ["10"]
        + emptvec
        + ["100"]
        + emptvec
    )
    return ticks, ticklabels


# function to find the stem (longest
# common substring) from the string array


def findstem(arr):
    # Determine size of the array
    n = len(arr)
    # Take first word from array
    # as reference
    s = arr[0]
    l = len(s)
    res = ""
    for i in range(l):
        for j in range(i + 2, l + 1):  # lenth at least 2
            # generating all possible substrings
            # of our reference string arr[0] i.e s
            stem = s[i:j]
            k = 1
            for k in range(n):
                # Check if the generated stem is
                # common to all words
                if stem not in arr[k]:
                    break
            # If current substring is present in
            # all strings and its length is greater
            # than current result
                if k + 1 == n and len(res) < len(stem):
                    res = stem
    return res


def findregexp(fnames):
    # use first name as a template
    baseStr = [fnames[0]]
    arr = baseStr + fnames

    # Function call
    stem = findstem(arr)
    stems = []
    while stem:
        # keep finding stems and replacing them in the template with stars. Can probably add an *ignore stars* to the stem finder
        baseStr = [baseStr[0].replace(stem, "*")]
        arr = baseStr + [s.replace(stem, "") for s in arr[1:]]
        # arr = baseStr+fnames
        stems.append(stem)
        stem = findstem(arr)

    # make a list of stars as lont as the OG template
    stars = ["*"] * len(fnames[0])
    # replace all stems in the list
    for s in stems:
        stars[fnames[0].find(s) : fnames[0].find(s) + len(s)] = s
    # remove repeating *s
    superStars = []
    superStars.append(stars[0])
    for s in stars[1:]:
        if not s == "*":
            superStars.append(s)
        elif not superStars[-1] == s:
            superStars.append(s)
    globExp = "".join(superStars)

    return globExp


def extractFieldsByRegex(globExp, fnames):
    import re

    matches = []
    for f in fnames:
        match = re.findall(
            globExp.replace("\\", "/").replace("*", "(.*)"), f.replace("\\", "/")
        )
        if not len(match) == 1:
            print("Non unique matches from regexp")
            break
        matches.append(match[0])
    return matches



def is_pathname_valid(pathname: str) -> bool:
    '''
    `True` if the passed pathname is a valid pathname for the current OS;
    `False` otherwise.
    '''
    import errno, os, sys
    ERROR_INVALID_NAME = 123
    try:
        if not isinstance(pathname, str) or not pathname:
            return False
        _, pathname = os.path.splitdrive(pathname)
        root_dirname = os.environ.get('HOMEDRIVE', 'C:') \
            if sys.platform == 'win32' else os.path.sep
        assert os.path.isdir(root_dirname)
        root_dirname = root_dirname.rstrip(os.path.sep) + os.path.sep
        for pathname_part in pathname.split(os.path.sep):
            try:
                os.lstat(root_dirname + pathname_part)
            except OSError as exc:
                if hasattr(exc, 'winerror'):
                    if exc.winerror == ERROR_INVALID_NAME:
                        return False
                elif exc.errno in {errno.ENAMETOOLONG, errno.ERANGE}:
                    return False
    except TypeError as exc:
        return False
    else:
        return True

def is_path_creatable(pathname: str) -> bool:
    '''
    `True` if the current user has sufficient permissions to create the passed
    pathname; `False` otherwise.
    '''
    import os
    # Parent directory of the passed path. If empty, we substitute the current
    # working directory (CWD) instead.
    dirname = os.path.dirname(pathname) or os.getcwd()
    return os.access(dirname, os.W_OK)

def is_path_exists_or_creatable(pathname: str) -> bool:
    '''
    `True` if the passed pathname is a valid pathname for the current OS _and_
    either currently exists or is hypothetically creatable; `False` otherwise.
    This function is guaranteed to _never_ raise exceptions.
    '''
    import os
    try:
        # To prevent "os" module calls from raising undesirable exceptions on
        # invalid pathnames, is_pathname_valid() is explicitly called first.
        return is_pathname_valid(pathname) and (
            os.path.exists(pathname) or is_path_creatable(pathname))
    # Report failure on non-fatal filesystem complaints (e.g., connection
    # timeouts, permissions issues) implying this path to be inaccessible. All
    # other exceptions are unrelated fatal issues and should not be caught here.
    except OSError:
        return False