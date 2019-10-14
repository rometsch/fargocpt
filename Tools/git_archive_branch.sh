#!/usr/bin/env bash
if [[ -z "$(git branch | grep $1)" ]]; then
	echo "Branch $1 does not exist!"
	exit 1
fi
git tag archive/$1 $1
git branch -D $1
git branch -d -r origin/$1
git push --tags
git push origin :$1
