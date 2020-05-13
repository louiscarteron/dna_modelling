import { processInput } from "../scripts/utils";

onmessage = function(e) {
  console.log('Worker: Message received from main script');
  console.log(e.data);
  console.log(processInput);
  random();
  postMessage("hello");
}

const random = () => {
  console.log("In the random function");
}
