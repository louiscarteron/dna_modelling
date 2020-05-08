import React from "react";
import { makeStyles } from "@material-ui/core/styles";
import Container from "@material-ui/core/Container";
import FileUpload from "../components/FileUpload";
import Summary from "../components/Summary";

const useStyles = makeStyles((theme) => ({
  mainGrid: {
    marginTop: theme.spacing(3)
  },
  mainContent: {
    backgroundColor: theme.palette.background.content,
    padding: theme.spacing(8),
    paddingTop: theme.spacing(4),
    paddingBottom: theme.spacing(4),
    minHeight: "calc(100vh - 64px)"
  }
}))

const Main = () => {

  const classes = useStyles();

  return (
    <Container maxWidth="lg" className={classes.mainContent}>
      <Summary/> 
      <FileUpload />
    </Container>
  )

}

export default Main;