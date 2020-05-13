import React, { useCallback, useState } from "react";
import { useDropzone } from "react-dropzone";
import Typography from "@material-ui/core/Typography";

import TextField from '@material-ui/core/TextField';
import { makeStyles } from "@material-ui/core/styles";
import Button from "@material-ui/core/Button";
import Card from '@material-ui/core/Card';
import CardActions from '@material-ui/core/CardActions';
import CardContent from '@material-ui/core/CardContent';

import ModelWorker from '../scripts/modelling.worker.js';

const useStyles = makeStyles(theme => ({
  textWrapper: {
    margin: theme.spacing(1),
    width: 400,
    margin: "auto",
    display: "flex",
    flexDirection: "column"
  },
  textField: {
    width: "100%"
  },
  runButton: {
    marginLeft: "auto",
    marginRight: theme.spacing(8)
  },
  uploadButton: {
    marginLeft: theme.spacing(8)
  },
  card: {
    margin: "auto",
    minWidth: 275,
    maxWidth: "50vw",
    height: '100%',
    display: 'flex',
    flexDirection: 'column',
    backgroundColor: theme.palette.background.paper,
    border: "none"
  },
  hiddenInput: {
    display: "none"
  }
}));

const FileUpload = (props) => {

  const classes = useStyles();

  const [inputs, setInputs] = useState("");

  const onDrop = useCallback((acceptedFiles) => {
    acceptedFiles.forEach((file) => {
      const reader = new FileReader();

      reader.onabort = () => console.log("Aborted");
      reader.onerror = () => console.log("Error");

      // TODO: Actually do something with the data. 
      reader.onload = () => {
        const bnyStr = reader.result;
        setInputs(bnyStr);
      }

      //reader.readAsArrayBuffer(file);
      reader.readAsBinaryString(file);
    });

  }, []);

  const { getRootProps, getInputProps } = useDropzone({onDrop: onDrop, noDrag: true});

  const runSimulation = () => {
    if (window.Worker) {
      const worker = new ModelWorker();
      worker.postMessage(inputs);

      worker.onmessage = (e) => {
        console.log("Recieved message");
        console.log(e.data);
      }

    } else {
      //Todo: Run blocking script instead of non-blocking
      console.log("Your browser doesn't support web workers.");
    }
  }

  return (
    <Card className={classes.card} elevation={0}>
      <CardContent>
        <Typography variant="body1">
          Upload a file and run a simulation on the inputs
        </Typography>
        <TextField
          className={classes.textField}
          rowsMax={20}
          placeholder="Drag file here or upload"
          multiline
          value={inputs}
          variant="outlined"
        /> 
      </CardContent>
      <CardActions disableSpacing>
        
        <div {...getRootProps()}>
          <input {...getInputProps()} className={classes.hiddenInput} id="contained-button-file"/>
          <label htmlFor="contained-button-file">
            <Button color="secondary" variant="outlined" className={classes.uploadButton}>
              Upload File
            </Button>
          </label>
        </div>

        <Button color="secondary" variant="outlined" className={classes.runButton} onClick={runSimulation}>
          Run simulation
        </Button>

      </CardActions>
    </Card>
  );
}

export default FileUpload;