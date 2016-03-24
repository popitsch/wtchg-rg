#!/bin/bash
export MAVEN_OPTS=-Xss2m
mvn $1 clean assembly:assembly -DskipTests=true -Denv=dev -U
mkdir -p bin
cp target/wtchg-rg-1.0-jar-with-dependencies.jar bin/wtchg-rg-1.0.jar
