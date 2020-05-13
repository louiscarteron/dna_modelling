export default () => {

  onmessage = function(e) {
    console.log('Worker: Message received from main script');
    console.log(e.data);
    random();
    postMessage("hello");
  }

  const random = () => {
    console.log("In the random function");
  }

}