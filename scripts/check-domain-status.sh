#!/bin/bash

echo "🌐 Domain Setup Status for queb.space"
echo "======================================"
echo ""

# Check DNS Zone Status
echo "📍 DNS Zone Status:"
if gcloud dns managed-zones describe queb-space-zone --format="value(name,dnsName,description)" 2>/dev/null; then
    echo "✅ DNS zone exists"
else
    echo "❌ DNS zone not found"
fi
echo ""

# Check DNS Records
echo "📋 DNS Records:"
if gcloud dns record-sets list --zone=queb-space-zone --format="table(name,type,ttl,rrdatas)" 2>/dev/null; then
    echo "✅ DNS records configured"
else
    echo "❌ Cannot list DNS records"
fi
echo ""

# Check Domain Verification Status
echo "🔍 Domain Verification Status:"
VERIFIED_DOMAIN=$(gcloud domains list-user-verified --filter="domain:queb.space" --format="value(domain)" 2>/dev/null)
if [ -n "$VERIFIED_DOMAIN" ]; then
    echo "✅ Domain verified: $VERIFIED_DOMAIN"
else
    echo "⏳ Domain verification pending"
    echo "   → Go to: https://search.google.com/search-console/welcome"
    echo "   → Add queb.space and verify ownership"
fi
echo ""

# Check Domain Mappings
echo "🔗 Domain Mappings:"
MAPPINGS=$(gcloud beta run domain-mappings list --region=us-central1 --format="table(name,domain,status)" 2>/dev/null)
if [ -n "$MAPPINGS" ]; then
    echo "$MAPPINGS"
else
    echo "⏳ No domain mappings found"
    echo "   → Run after domain verification:"
    echo "   → gcloud beta run domain-mappings create --service=molecular-analysis --domain=queb.space --region=us-central1"
    echo "   → gcloud beta run domain-mappings create --service=molecular-analysis --domain=www.queb.space --region=us-central1"
fi
echo ""

# Check Cloud Function Status
echo "☁️ Cloud Function Status:"
FUNCTION_STATUS=$(gcloud functions describe molecular-analysis --region=us-central1 --format="value(name,status,httpsTrigger.url)" 2>/dev/null)
if [ -n "$FUNCTION_STATUS" ]; then
    echo "✅ Cloud Function: $FUNCTION_STATUS"
else
    echo "❌ Cloud Function not found"
fi
echo ""

# Check if nameservers are properly configured
echo "🔧 Nameserver Configuration:"
echo "Current Google Cloud nameservers:"
gcloud dns managed-zones describe queb-space-zone --format="value(nameServers)" 2>/dev/null | tr ';' '\n' | sed 's/^/   /'
echo ""

# Check DNS propagation
echo "🌍 DNS Propagation Check:"
if command -v dig &> /dev/null; then
    echo "Checking queb.space nameservers:"
    dig queb.space NS +short 2>/dev/null | sed 's/^/   /' || echo "   ❌ Cannot check nameservers"
    
    echo "Checking www.queb.space:"
    dig www.queb.space A +short 2>/dev/null | sed 's/^/   /' || echo "   ❌ Cannot resolve www.queb.space"
else
    echo "⚠️ dig command not available - install bind-utils to check DNS propagation"
fi
echo ""

# Provide actionable next steps
echo "✅ Next Steps:"
echo "1. Update nameservers at Namecheap to Google Cloud nameservers above"
echo "2. Wait 24-48 hours for DNS propagation"
echo "3. Complete domain verification in Google Search Console"
echo "4. Run domain mapping commands shown above"
echo "5. Test: https://queb.space and https://www.queb.space"
echo ""

# Check if we can reach the current function URL
echo "🔗 Current Function URL Test:"
FUNCTION_URL=$(gcloud functions describe molecular-analysis --region=us-central1 --format="value(httpsTrigger.url)" 2>/dev/null)
if [ -n "$FUNCTION_URL" ]; then
    echo "Testing: $FUNCTION_URL"
    if curl -s --max-time 5 "$FUNCTION_URL" > /dev/null 2>&1; then
        echo "✅ Function is responding"
    else
        echo "❌ Function not responding (may be cold start)"
    fi
else
    echo "⚠️ Cannot get function URL - checking function status..."
    gcloud functions describe molecular-analysis --region=us-central1 --format="value(status)" 2>/dev/null || echo "❌ Function not found"
fi 