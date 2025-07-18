#!/bin/bash
# Development script - run tests and start server with nodemon
jest --testPathPattern=unit.test.js --verbose --silent && nodemon server.js 