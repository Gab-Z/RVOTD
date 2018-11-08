/*
global.eval = function () {
  throw new Error(`Sorry, this app does not support window.eval().`)
}*/

const electron = require( 'electron' );
const { app, BrowserWindow, ipcMain, session } = require( 'electron' );

const path = require( 'path' );
const url = require( 'url' );
let mainWindow;

function createWindow () {

  // Create the browser window.
  mainWindow = new BrowserWindow( {
    width: 1750,
    height: 1100

  } );
  // and load the index.html of the app.
  mainWindow.webContents.openDevTools();

  mainWindow.loadURL( url.format( {
    pathname: path.join( __dirname, 'html/index.html' ),
    protocol: 'file:',
    slashes: true
  }))

  mainWindow.on( 'closed', () => {
    mainWindow = null;
  })

}
app.on( 'ready', createWindow );

app.on( 'window-all-closed', () => {
  if ( process.platform !== 'darwin' ) {
    app.quit();
  }
});

app.on( 'activate', () => {
  if (mainWindow === null) {
    createWindow();
  }
});
