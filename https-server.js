const https = require("https");
const fs = require("fs");
const path = require("path");
const os = require("os");
const { execSync } = require("child_process");

class HttpsServer {
  constructor(app, port = 3001) {
    this.app = app;
    this.port = port;
    this.localIP = this.getLocalIPAddress();
  }

  getLocalIPAddress() {
    const interfaces = os.networkInterfaces();
    for (const name of Object.keys(interfaces)) {
      for (const iface of interfaces[name]) {
        if (iface.family === "IPv4" && !iface.internal) {
          return iface.address;
        }
      }
    }
    return "127.0.0.1";
  }

  generateSelfSignedCert() {
    const certDir = path.join(__dirname, "certs");
    const keyPath = path.join(certDir, "key.pem");
    const certPath = path.join(certDir, "cert.pem");
    const configPath = path.join(certDir, "openssl.conf");

    if (fs.existsSync(keyPath) && fs.existsSync(certPath)) {
      console.log("✅ Using existing SSL certificates");
      return { key: fs.readFileSync(keyPath), cert: fs.readFileSync(certPath) };
    }

    if (!fs.existsSync(certDir)) fs.mkdirSync(certDir, { recursive: true });

    console.log(
      "🔐 Generating self-signed SSL certificates for development...",
    );

    try {
      const configContent = `[req]
distinguished_name = req_distinguished_name
req_extensions = v3_req
prompt = no

[req_distinguished_name]
C = US
ST = Development
L = Local
O = MolecularReality
CN = localhost

[v3_req]
basicConstraints = CA:FALSE
keyUsage = nonRepudiation, digitalSignature, keyEncipherment
subjectAltName = @alt_names

[alt_names]
DNS.1 = localhost
DNS.2 = *.localhost
IP.1 = 127.0.0.1
IP.2 = 0.0.0.0
IP.3 = ${this.localIP}
IP.4 = ::1
`;

      fs.writeFileSync(configPath, configContent);

      execSync(`openssl genrsa -out ${keyPath} 2048`, { stdio: "inherit" });
      execSync(
        `openssl req -new -x509 -key ${keyPath} -out ${certPath} -days 365 -config ${configPath}`,
        { stdio: "inherit" },
      );

      fs.unlinkSync(configPath);

      console.log("✅ SSL certificates generated successfully!");
      return { key: fs.readFileSync(keyPath), cert: fs.readFileSync(certPath) };
    } catch (err) {
      console.error("❌ Failed to generate SSL certificates:", err.message);
      return null;
    }
  }

  start() {
    const credentials = this.generateSelfSignedCert();
    if (credentials) {
      const server = https.createServer(credentials, this.app);
      server.listen(this.port, "0.0.0.0", () => {
        console.log(
          `🔒 HTTPS server running on https://${this.localIP}:${this.port}`,
        );
      });
      return server;
    } else {
      console.log("⚠️ HTTPS not started — use ngrok if needed");
      return null;
    }
  }
}

module.exports = HttpsServer;
