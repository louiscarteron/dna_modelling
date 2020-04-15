import React from 'react';
import logo from './logo.svg';
import './App.css';
import 'normalize.css';

function App() {
  return (
    <div className="App">
      <header className="App-header">
        <img src={logo} className="App-logo" alt="logo" />
        <p>
          Edit <code>src/App.js</code> and save to reload.
          What if I don't want to huh?
        </p>
        <a
          className="App-link"
          href="https://reactjs.org"
          target="_blank"
          rel="noopener noreferrer"
        >
          Learn React. Or maybe not 😋😋
        </a>
      </header>
    </div>
  );
}

export default App;
