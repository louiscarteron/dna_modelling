import React, { useState } from 'react';
import logo from './logo.svg';
import './App.css';
import 'normalize.css';
import CssBaseline from "@material-ui/core/CssBaseline";
import Appbar from "./components/Appbar";

import { createMuiTheme, ThemeProvider } from "@material-ui/core";
import Main from './pages/Main';

// TODO: Find actual theme colours 
const lightTheme = createMuiTheme({
  palette: {
    primary: {
      main: '#111'
    },
    secondary: {
      main: '#ff2266'
    },
    text: {
      primary: "#000",
      secondary: "#000",
      disabled: '#ff0000'
    },
    background: {
      //default: "#bbb",
      //paper: "#f1f1f1",
      default: "#fff",
      paper: "#f0f0f0",
      content: "#f0f0f0"
    }
  }
});

const darkTheme = createMuiTheme({
  palette: {
    primary: {
      main: '#111'
    },
    secondary: {
      main: '#ff2266'
    },
    text: {
      primary: "#fff",
      secondary: "#fff",
      disabled: '#ff0000'
    },
    background: {
      //default: "#bbb",
      //paper: "#f1f1f1",
      default: "#111",
      paper: "#0f0f0f",
      content: "#0f0f0f"
    }
  }
});

function App() {

  const [theme, setDarkTheme] = useState(false);

  return (
    <ThemeProvider theme={theme ? darkTheme : lightTheme}>
      <CssBaseline />
      <Appbar switchTheme={setDarkTheme}/>
      <Main/>
    </ThemeProvider>
  );
}

export default App;
