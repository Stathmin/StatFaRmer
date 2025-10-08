# StatFaRmer Security & Deployment Guide

> **Status**: Planned project - comprehensive security and deployment documentation

## 🔒 Security Configuration

### Current Security Status

#### ✅ **Good Security Practices:**
- **Project Filtering**: `deploy_config.R` restricts visible projects in public deployment
- **Consistent Usage**: Configuration files are properly loaded and used
- **Data Isolation**: Only allowed projects are accessible

#### ⚠️ **Security Concerns:**
- **Network Binding**: Default `0.0.0.0` exposes app to all network interfaces
- **No Authentication**: No user authentication system
- **No HTTPS**: No SSL/TLS encryption
- **No Rate Limiting**: Vulnerable to abuse

### 🛡️ Security Improvements Implemented

#### 1. **Secure Default Configuration**
```r
# config/global_config.R - Updated
SHINY_HOST <- '127.0.0.1'  # Localhost only - secure default
```

#### 2. **Environment-Based Security**
```r
# config/security_config.R - New
# Production: 127.0.0.1 only, authentication enabled
# Development: Local networks allowed, no auth
```

#### 3. **Project Access Control**
```r
# config/deploy_config.R - Existing
allowed_projects <- c('project_NO3', 'project_soy_2024-05')
```

## 🚀 Deployment Recommendations

### **Development Environment**
```bash
# Safe for local development
export STATFARMER_ENV=development
Rscript launch_statfarmer.R wizard
```

### **Production Deployment**

#### **Option 1: Reverse Proxy (Recommended)**
```bash
# Use nginx/Apache as reverse proxy
# Bind Shiny to localhost only
export STATFARMER_ENV=production
Rscript launch_statfarmer.R app
```

**Nginx Configuration:**
```nginx
server {
    listen 443 ssl;
    server_name your-domain.com;
    
    ssl_certificate /path/to/cert.pem;
    ssl_certificate_key /path/to/key.pem;
    
    location / {
        proxy_pass http://127.0.0.1:3839;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }
}
```

#### **Option 2: Docker Deployment**
```dockerfile
FROM rocker/shiny:latest

# Copy application
COPY . /srv/shiny-server/statfarmer/

# Set environment
ENV STATFARMER_ENV=production

# Expose only to localhost
EXPOSE 3839

# Run with localhost binding
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/statfarmer', host='127.0.0.1', port=3839)"]
```

### **Security Checklist**

#### **Before Public Deployment:**
- [ ] Set `STATFARMER_ENV=production`
- [ ] Use reverse proxy with HTTPS
- [ ] Configure firewall (only allow proxy access)
- [ ] Review `allowed_projects` in `deploy_config.R`
- [ ] Enable access logging
- [ ] Set up monitoring/alerting
- [ ] Regular security updates

#### **Network Security:**
- [ ] Bind to `127.0.0.1` only (not `0.0.0.0`)
- [ ] Use HTTPS with valid certificates
- [ ] Configure firewall rules
- [ ] Monitor access logs
- [ ] Implement rate limiting

#### **Data Security:**
- [ ] Restrict project access via `deploy_config.R`
- [ ] Regular backups of data
- [ ] Secure file permissions
- [ ] No sensitive data in public projects

## 🔧 Configuration Files

### **config/global_config.R**
- System-wide settings
- Secure defaults (localhost binding)
- Processing parameters

### **config/deploy_config.R**
- Public deployment restrictions
- Project access control
- Data filtering

### **config/security_config.R** (New)
- Environment-based security
- Production vs development settings
- Network access controls

## 📊 Current Status

| Security Aspect | Status | Notes |
|----------------|--------|-------|
| Project Filtering | ✅ Good | `deploy_config.R` working |
| Network Binding | ✅ Fixed | Changed to localhost |
| Authentication | ⚠️ Missing | Consider adding |
| HTTPS | ⚠️ Missing | Use reverse proxy |
| Rate Limiting | ⚠️ Missing | Consider adding |
| Access Logging | ⚠️ Missing | Consider adding |

## 🎯 Next Steps

1. **Immediate**: Use localhost binding (already implemented)
2. **Short-term**: Set up reverse proxy with HTTPS
3. **Medium-term**: Add authentication system
4. **Long-term**: Implement comprehensive security monitoring

---

**Remember**: Security is an ongoing process. Regularly review and update your security configuration!
