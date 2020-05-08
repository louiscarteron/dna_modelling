import React from "react";

import Typography from "@material-ui/core/Typography";

import { makeStyles } from "@material-ui/core/styles";
import Card from '@material-ui/core/Card';
import CardContent from '@material-ui/core/CardContent';

const useStyles = makeStyles(theme => ({
  card: {
    margin: "auto",
    minWidth: 275,
    maxWidth: "50vw",
    height: '100%',
    display: 'flex',
    flexDirection: 'column',
    backgroundColor: theme.palette.background.paper,
    border: "none"
  }
}));

const Summary = () => {

  const classes = useStyles();

  return(

    <Card className={classes.card} elevation={0}>
      <CardContent>
        <Typography variant="body1">
          DNA storage blah blah blah
          Lorem ipsum dolor, sit amet consectetur adipisicing elit. Ratione, harum fugiat sint quam aspernatur facilis ea incidunt aliquam facere hic nesciunt? Reiciendis doloremque nemo, adipisci sit iusto omnis inventore ipsum!
        </Typography>
      </CardContent>
    </Card>

  );

}

export default Summary;