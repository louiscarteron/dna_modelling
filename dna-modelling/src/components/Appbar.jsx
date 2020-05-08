import React, { useState } from "react";
import AppBar from "@material-ui/core/AppBar";
import Typography from "@material-ui/core/Typography";
import { withStyles } from "@material-ui/core/styles";
import Toolbar from "@material-ui/core/Toolbar";

import ToggleButton from '@material-ui/lab/ToggleButton';
import Brightness7Icon from '@material-ui/icons/Brightness7';
import Brightness4Icon from '@material-ui/icons/Brightness4';

const styles = {
  root: {
    flexGrow: 1,
  },
  grow: {
    flexGrow: 1
  }
}

// Double Appbar to take into account width of appbar and offset page content accordingly 
const Appbar = (props) => {

  const { classes, switchTheme } = props;

  const updateTheme = () => {
    switchTheme(!selected);
    setSelected(!selected);
  }

  const [selected, setSelected] = useState(false);

  return (
    <div className={classes.root}>
      <AppBar color="secondary">
        <Toolbar>
          <Typography variant="h6" color="inherit" className={classes.grow}>
            DNA Modelling
          </Typography>
          <ToggleButton
            value="check"
            selected={selected}
            onChange={updateTheme}
          >
            {selected ? <Brightness7Icon /> : <Brightness4Icon/>}
          </ToggleButton>
        </Toolbar>
      </AppBar>
      <Toolbar/>
    </div>
  );
}

export default withStyles(styles)(Appbar);