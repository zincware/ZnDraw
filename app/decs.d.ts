declare module "@json-editor/json-editor" {
  const JSONEditor: any;
  export default JSONEditor;
}
declare module "znsocket" {
  function createClient(...args: any[]): any;
  const Client: any;
  const List: any;
  const Dict: any;
  export { createClient, List, Dict, Client };
}
